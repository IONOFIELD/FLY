"""
Connectome-constrained LIF network in Brian2.

Effective chemical kick per synapse = w_syn_mV * w_scale = 0.275 * 0.3 = 0.0825 mV.
0.275 mV is Shiu et al. 2024's unitary weight on FlyWire; 0.3 is the MaleCNS
rescale found by sweep (VNC + hubs make the raw value unstable).

Neuron:  dv/dt = (V_rest - v + g + I_gap - a)/tau_m ; dg/dt = -g/tau_syn ;
         da/dt = -a/tau_adapt ; on spike a += b_adapt (0 for all but overridden types)
Chemical synapse: g_post += sign * n_synapses * w_syn * w_scale
Electrical synapse: I_gap_post = g_gap*(v_pre - v_post)  (summed) plus a
         spikelet v_post += A_spikelet with the measured 0.8 ms latency.
Parameters from Shiu et al. 2024 unless stated; w_scale fitted on MaleCNS.
"""
from dataclasses import dataclass, field
import numpy as np
import pandas as pd
from brian2 import (NeuronGroup, PoissonGroup, Synapses, SpikeMonitor, StateMonitor, Network,
                    ms, mV, Hz, second, defaultclock, prefs, seed as b2seed)

from .graph import electrical_pairs, drop_mixed_chemical, pre_class, INTRINSIC_OVERRIDES


@dataclass
class LIFParams:
    v_rest_mV: float = -52.0
    v_thresh_mV: float = -45.0
    v_reset_mV: float = -52.0
    tau_m_ms: float = 20.0
    tau_syn_ms: float = 5.0
    refractory_ms: float = 2.2
    w_syn_mV: float = 0.275          # Shiu et al. 2024 per-synapse kick
    w_scale: float = 0.3             # fitted on MaleCNS (sweep.py)
    tau_adapt_ms: float = 100.0
    gap_delay_ms: float = 0.8        # Tanouye & Wyman 1980
    chem_delay_ms: float = 1.8       # Shiu et al. 2024
    dt_ms: float = 0.1
    seed: int = 0
    # class-wise gain on top of w_scale, keyed by presynaptic class (graph.pre_class)
    class_gains: dict = field(default_factory=lambda: {"sensory": 1.0, "relay": 1.0, "local": 1.0})


class CNSModel:
    def __init__(self, neurons, edges, stim_types, params=LIFParams(),
                 electrical=True, codegen="cython", monitor_types=()):
        prefs.codegen.target = codegen
        defaultclock.dt = params.dt_ms * ms
        b2seed(params.seed)
        self.p = params
        self.neurons = neurons
        self.meta = neurons.set_index("bodyId")

        if electrical:
            edges = drop_mixed_chemical(edges, neurons)
        stim = neurons[neurons["type"].isin(stim_types)][["bodyId", "type"]].reset_index(drop=True)
        self.stim = stim
        self.stim_index = pd.Series(np.arange(len(stim)), index=stim["bodyId"].values)
        all_ids = pd.unique(pd.concat([edges["pre"], edges["post"]]))
        self.lif_ids = np.array([b for b in all_ids if b not in self.stim_index.index])
        self.lif_index = pd.Series(np.arange(len(self.lif_ids)), index=self.lif_ids)

        p = params
        eqs = """
        dv/dt = (v_rest - v + g + I_gap - a) / tau_m : volt (unless refractory)
        dg/dt = -g / tau_syn : volt
        da/dt = -a / tau_adapt : volt
        I_gap : volt
        b_adapt : volt
        """
        ns = dict(v_rest=p.v_rest_mV * mV, tau_m=p.tau_m_ms * ms, tau_syn=p.tau_syn_ms * ms,
                  tau_adapt=p.tau_adapt_ms * ms, v_thresh=p.v_thresh_mV * mV,
                  v_reset=p.v_reset_mV * mV)
        G = NeuronGroup(len(self.lif_ids), eqs, threshold="v > v_thresh",
                        reset="v = v_reset; a += b_adapt",
                        refractory=p.refractory_ms * ms, method="euler", namespace=ns)
        G.v = p.v_rest_mV * mV
        G.b_adapt = 0 * mV
        self.overrides_applied = []
        for typ, par, val, src in INTRINSIC_OVERRIDES:
            idx = [self.lif_index[b] for b in self.meta.index[self.meta["type"] == typ]
                   if b in self.lif_index.index]
            if par == "b_adapt_mV" and idx:
                G.b_adapt[idx] = val * mV
                self.overrides_applied.append(dict(type=typ, parameter=par, value=val, source=src))

        w_unit = p.w_syn_mV * p.w_scale * mV
        gain_of = self.meta["superclass"].map(pre_class).map(p.class_gains).fillna(1.0)
        rec = edges[edges["pre"].isin(self.lif_index.index) & edges["post"].isin(self.lif_index.index)]
        S_rec = Synapses(G, G, "w : volt", on_pre="g_post += w")
        S_rec.connect(i=self.lif_index[rec["pre"]].values, j=self.lif_index[rec["post"]].values)
        S_rec.w = (rec["weight"] * rec["sign"] * gain_of.loc[rec["pre"]].values).values * w_unit
        S_rec.delay = p.chem_delay_ms * ms

        P = PoissonGroup(len(stim), rates=0 * Hz)
        se = edges[edges["pre"].isin(self.stim_index.index) & edges["post"].isin(self.lif_index.index)]
        objs = [G, P, S_rec]
        if len(se):
            S_stim = Synapses(P, G, "w : volt", on_pre="g_post += w")
            S_stim.connect(i=self.stim_index[se["pre"]].values, j=self.lif_index[se["post"]].values)
            S_stim.w = (se["weight"] * se["sign"] * gain_of.loc[se["pre"]].values).values * w_unit
            S_stim.delay = p.chem_delay_ms * ms
            objs.append(S_stim)
        else:
            print("[model] WARNING: stimulated types have no outgoing edges to LIF neurons")
        self.electrical_table = pd.DataFrame()
        if electrical:
            ep_all = electrical_pairs(neurons)
            ep_all = ep_all[ep_all.post.isin(self.lif_index.index)]
            # LIF -> LIF: ohmic coupling + spikelet
            ep = ep_all[ep_all.pre.isin(self.lif_index.index)]
            if len(ep):
                S_gap = Synapses(G, G,
                                 "I_gap_post = gg * (v_pre - v_post) : volt (summed)\n spk : volt\n gg : 1",
                                 on_pre="v_post += spk")
                S_gap.connect(i=self.lif_index[ep.pre].values, j=self.lif_index[ep.post].values)
                S_gap.spk = ep.spikelet_mV.values * mV
                S_gap.gg = ep.g_gap.values
                S_gap.delay = p.gap_delay_ms * ms
                objs.append(S_gap)
            # stimulus source -> LIF: spikelet only (sources have no membrane voltage)
            eps = ep_all[ep_all.pre.isin(self.stim_index.index)]
            if len(eps):
                S_gap_s = Synapses(P, G, "spk : volt", on_pre="v_post += spk")
                S_gap_s.connect(i=self.stim_index[eps.pre].values, j=self.lif_index[eps.post].values)
                S_gap_s.spk = eps.spikelet_mV.values * mV
                S_gap_s.delay = p.gap_delay_ms * ms
                objs.append(S_gap_s)
            self.electrical_table = ep_all

        self.spikes = SpikeMonitor(G)
        objs.append(self.spikes)
        self.vmon = None
        mon_idx = [self.lif_index[b] for t in monitor_types
                   for b in self.meta.index[self.meta["type"] == t] if b in self.lif_index.index]
        if mon_idx:
            self.vmon = StateMonitor(G, "v", record=mon_idx)
            self.vmon_ids = self.lif_ids[mon_idx]
            objs.append(self.vmon)
        self.G, self.P = G, P
        self.net = Network(*objs)
        self.n_chem = len(rec)
        self.n_stim_syn = len(se)

    def store(self):
        self.net.store("base")

    def restore(self):
        self.net.restore("base")

    def readout_index(self, typ):
        return [self.lif_index[b] for b in self.meta.index[self.meta["type"] == typ]
                if b in self.lif_index.index]

    def run(self, duration_ms):
        self.net.run(duration_ms * ms)

    def set_stim_rates(self, rates_hz):
        self.P.rates = np.asarray(rates_hz) * Hz

    @property
    def t_ms(self):
        return float(self.net.t / ms)

    def spike_frame(self):
        t = np.array(self.spikes.t / ms)
        i = np.array(self.spikes.i)
        df = pd.DataFrame({"t_ms": t, "bodyId": self.lif_ids[i]})
        df["type"] = df["bodyId"].map(self.meta["type"])
        df["superclass"] = df["bodyId"].map(self.meta["superclass"])
        return df

    def peak_depolarization_mV(self, onsets, window_ms):
        """Max membrane deflection above rest (mV) of monitored cells in each window."""
        if self.vmon is None:
            return np.array([])
        t = np.array(self.vmon.t / ms); v = np.array(self.vmon.v / mV)
        out = []
        for on in onsets:
            m = (t >= on) & (t < on + window_ms)
            out.append((v[:, m].max(axis=1) - self.p.v_rest_mV).mean() if m.any() else np.nan)
        return np.array(out)

    def population_rate_hz(self):
        return len(self.spikes.i) / len(self.lif_ids) / (self.t_ms / 1000.0)
