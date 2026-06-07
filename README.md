# cech-complex

> **Čech complex from point cloud via ball intersection nerve**

[![crates.io](https://img.shields.io/crates/v/cech-complex.svg)](https://crates.io/crates/cech-complex)
[![license](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)

The Čech complex is the nerve of balls of radius ε around each point. Unlike the Vietoris-Rips complex (which is an approximation), the Čech complex has a homotopy equivalence guarantee via the Nerve Theorem.

## Construction

For radius ε, the Čech complex Cech(ε) has:
- A vertex for each point
- A k-simplex {p₀, ..., pₖ} if ∩ᵢ B(pᵢ, ε) ≠ ∅

## Advantages over VR

- **Homotopy equivalence**: Nerve theorem guarantees topological correctness
- **Smaller**: Fewer simplices than VR at the same scale
- **Exact**: No approximation — captures true topology

## Tradeoff

The Čech complex is harder to compute (requires checking all intersections), but is topologically exact.

## Installation

```toml
[dependencies]
cech-complex = "0.1.0"
```

## License

MIT © [SuperInstance](https://github.com/SuperInstance)

---

*Part of the [Exocortex](https://github.com/SuperInstance/exocortex) project — persistent cognitive substrate for multi-agent systems.*
