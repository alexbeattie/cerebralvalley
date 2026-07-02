"""A compact sequence -> activity CNN (Basset/DeepSTARR-like), CPU-trainable.

Why this shape: convolutions over one-hot DNA are the standard MPRA-model inductive bias --
each filter is a learnable motif detector, and stacking a couple of them lets the head
compose motif presence into activity. Global average+max pooling makes the head
length-agnostic (so ISM can run on any window) and, with a planted-motif task, keeps the
model small enough to train in seconds on a laptop while still recovering the signal.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass

import torch
from torch import nn


@dataclass
class ModelConfig:
    """Architecture knobs, serialized into the model sidecar for reproducibility."""

    seq_length: int
    conv_channels: tuple[int, ...] = (64, 128)
    kernel_size: int = 9
    pool_size: int = 2
    mlp_hidden: int = 64
    dropout: float = 0.2
    in_channels: int = 4  # ACGT


class ActivityCNN(nn.Module):
    """One-hot (B, 4, L) -> scalar activity (B,).

    Conv blocks (Conv1d -> ReLU -> MaxPool) act as motif detectors; we then pool over the
    length axis with both mean and max (mean = overall density, max = strongest hit) and
    map the concatenated summary to a scalar with a small MLP.
    """

    def __init__(self, config: ModelConfig) -> None:
        super().__init__()
        self.config = config

        blocks: list[nn.Module] = []
        in_ch = config.in_channels
        for out_ch in config.conv_channels:
            blocks.append(
                nn.Conv1d(in_ch, out_ch, config.kernel_size, padding=config.kernel_size // 2)
            )
            blocks.append(nn.ReLU())
            blocks.append(nn.MaxPool1d(config.pool_size))
            in_ch = out_ch
        self.conv = nn.Sequential(*blocks)

        # mean + max over length -> 2 * last channel count.
        feat = 2 * config.conv_channels[-1]
        self.head = nn.Sequential(
            nn.Linear(feat, config.mlp_hidden),
            nn.ReLU(),
            nn.Dropout(config.dropout),
            nn.Linear(config.mlp_hidden, 1),
        )

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        h = self.conv(x)  # (B, C, L')
        pooled = torch.cat([h.mean(dim=2), h.amax(dim=2)], dim=1)  # (B, 2C)
        return self.head(pooled).squeeze(-1)  # (B,)


def build_model(config: ModelConfig) -> ActivityCNN:
    return ActivityCNN(config)


def config_to_dict(config: ModelConfig) -> dict:
    d = asdict(config)
    d["conv_channels"] = list(config.conv_channels)  # JSON-friendly
    return d


def config_from_dict(d: dict) -> ModelConfig:
    d = dict(d)
    if "conv_channels" in d:
        d["conv_channels"] = tuple(d["conv_channels"])
    return ModelConfig(**d)
