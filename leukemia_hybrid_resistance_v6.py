
from dataclasses import dataclass, asdict
from typing import Tuple, Dict, Any, List
import numpy as np
from scipy.ndimage import convolve

@dataclass
class Params:
    grid_size: int = 100
    carrying_capacity_per_site: int = 2
    p_stem_div_base: float = 0.05
    p_diff_base: float = 0.08
    support_radius: int = 2
    support_source_rate: float = 0.35
    support_decay: float = 0.05
    support_diffusivity: float = 0.15
    chemo_kill_Ls: float = 0.65
    chemo_kill_Lr: float = 0.08
    chemo_kill_H: float = 0.12
    mutation_rate: float = 0.003
    max_days: int = 40
    seed: int = 42

class HybridMarrow:
    def __init__(self, p: Params):
        self.p = p
        self.rng = np.random.default_rng(p.seed)
        self.H = np.zeros((p.grid_size, p.grid_size), dtype=np.int32)
        self.Ls = np.zeros((p.grid_size, p.grid_size), dtype=np.int32)
        self.Lr = np.zeros((p.grid_size, p.grid_size), dtype=np.int32)
        self.S = np.zeros((p.grid_size, p.grid_size), dtype=np.float64)
        
        # مقداردهی اولیه سلول‌های نرمال و لوسمی
        mid = p.grid_size // 2
        self.H[mid-5:mid+5, mid-5:mid+5] = 1
        self.Ls[mid-2:mid+2, mid-2:mid+2] = 1

    def step_support(self):
        kernel = np.array([[0, 1, 0], [1, -4, 1], [0, 1, 0]], dtype=np.float64)
        lap = convolve(self.S, kernel, mode='wrap')
        prod = (self.H > 0).astype(np.float64) * self.p.support_source_rate
        self.S += self.p.support_diffusivity * lap + prod - self.p.support_decay * self.S
        np.clip(self.S, 0.0, 1.0, out=self.S)

    def run(self, snapshot_times: Tuple[int, ...] = (0,)):
        p = self.p
        H_hist, Ls_hist, Lr_hist = [], [], []
        total_leu_hist, total_cells_hist = [], []
        
        for day in range(p.max_days + 1):
            tot_h = int(self.H.sum())
            tot_ls = int(self.Ls.sum())
            tot_lr = int(self.Lr.sum())
            tot_leu = tot_ls + tot_lr
            tot_all = tot_h + tot_leu
            
            H_hist.append(tot_h)
            Ls_hist.append(tot_ls)
            Lr_hist.append(tot_lr)
            total_leu_hist.append(tot_leu)
            total_cells_hist.append(tot_all)
            
            if day == p.max_days:
                break
                
            self.step_support()
            
            # تغییرات تجمعی بر اساس احتمالات
            div_mask = (self.H + self.Ls + self.Lr) < (p.carrying_capacity_per_site * 4)
            if np.any(div_mask):
                h_grow = (self.H > 0) & (self.rng.random(self.H.shape) < p.p_stem_div_base)
                self.H[h_grow] += 1
                
                ls_grow = (self.Ls > 0) & (self.rng.random(self.Ls.shape) < p.p_stem_div_base * 1.2)
                self.Ls[ls_grow] += 1
                
                mut = (self.Ls > 0) & (self.rng.random(self.Ls.shape) < p.mutation_rate)
                self.Ls[mut] -= 1
                self.Lr[mut] += 1
                
        # تولید خروجی آبجکت نتیجه
        class Result:
            pass
        res = Result()
        res.H = H_hist
        res.Ls = Ls_hist
        res.Lr = Lr_hist
        res.total_leukemia = total_leu_hist
        res.total_cells = total_cells_hist
        res.final_positions = np.argwhere((self.H + self.Ls + self.Lr) > 0)
        res.final_cell_types = self.H * 1 + self.Ls * 2 + self.Lr * 3
        res.final_S = self.S
        return res
