# cech-complex

> **Čech complex from point clouds — exact topology via the Nerve Theorem**

[![crates.io](https://img.shields.io/crates/v/cech-complex.svg)](https://crates.io/crates/cech-complex)
[![docs.rs](https://docs.rs/cech-complex/badge.svg)](https://docs.rs/cech-complex)
[![license](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)

## What is the Čech Complex?

Given a set of points P in ℝᵈ and a radius r, place a ball of radius r around each point. The **Čech complex** is a simplicial complex where a k-simplex {p₀, ..., pₖ} exists if and only if the intersection of all k+1 balls is non-empty. In other words, a simplex exists when there exists at least one point in space that is within distance r of all vertices.

This is the **nerve** of the ball cover — the combinatorial record of which balls overlap. The critical mathematical property is the **Nerve Theorem**: when the underlying sets are convex (as balls in Euclidean space always are), the nerve is homotopy equivalent to the union of the sets. This means the Čech complex captures the **exact topology** of the point cloud at scale r — no approximation.

## Why Does This Matter?

The Čech complex is the "gold standard" for topological data analysis:

- **Topological exactness**: Unlike the Vietoris-Rips complex (which only approximates), the Čech complex is guaranteed to have the correct homotopy type by the Nerve Theorem
- **Smaller complexes**: At the same radius, the Čech complex has fewer simplices than the VR complex because it requires all-way intersection, not just pairwise
- **Filtration structure**: Sweeping r from 0 to ∞ produces a filtration — the basis for persistent homology computation
- **Theoretical foundation**: The Čech complex is the theoretical benchmark against which other complexes (VR, witness, alpha) are measured

Real-world applications:
- **Sensor networks**: Determine coverage holes — the Čech complex at the sensing radius exactly captures which areas are covered
- **Molecular biology**: Model protein structures where atoms are balls and the Čech complex captures the void/tunnel structure
- **Shape reconstruction**: Recover the topology of a surface from a point sample
- **Robotics**: Configuration space analysis — determine if a robot can navigate through obstacles

## Architecture

```
┌──────────────────────────────────────────────────────────────┐
│                   Čech Complex Pipeline                       │
│                                                              │
│  Point Cloud      Balls of radius r     Ball Intersections   │
│  ┌─────┐         ┌─────────────┐       ┌─────────────┐      │
│  │ p₀  │         │  ○       ○  │       │ p₀∩p₁ ≠ ∅  │      │
│  │ p₁  │ ────▶  │    ○   ○    │ ───▶  │ p₁∩p₂ ≠ ∅  │      │
│  │ p₂  │         │  ○       ○  │       │ p₀∩p₂ ≠ ∅  │      │
│  └─────┘         └─────────────┘       │ p₀∩p₁∩p₂?  │      │
│                                       └──────┬──────┘      │
│                                              │               │
│                    Nerve Construction         ▼               │
│                    ┌──────────────────────────────┐          │
│                    │      Čech Complex C(r)       │          │
│                    │                              │          │
│                    │   • p₀ ──── • p₁            │          │
│                    │    ╲       ╱  │              │          │
│                    │     ╲   ╱    │              │          │
│                    │      • p₂    │              │          │
│                    │                              │          │
│                    │  △ = {p₀,p₁,p₂} if triple   │          │
│                    │      intersection non-empty   │          │
│                    └──────────────────────────────┘          │
│                                                              │
│  ┌────────────────────────────────────────────────────────┐  │
│  │  Radius Sweep: r₁ < r₂ < r₃ → filtration C(r₁) ⊂ ... │  │
│  │  Connection radius: minimum r where all points connect │  │
│  └────────────────────────────────────────────────────────┘  │
└──────────────────────────────────────────────────────────────┘
```

## Quick Start

```rust
use cech_complex::{CechComplex, RadiusSweep, HomotopyEquivalence};

// Three points forming a triangle
let points = vec![
    vec![0.0, 0.0],
    vec![1.0, 0.0],
    vec![0.5, 0.866],  // equilateral triangle
];

// Build the Čech complex at radius 1.0, up to dimension 2
let complex = CechComplex::build(&points, 1.0, 2);
println!("Vertices: {}", complex.num_vertices());
println!("Edges: {}", complex.simplices_of_dim(1).len());
println!("Triangles: {}", complex.simplices_of_dim(2).len());

// Check connectivity
println!("Connected: {}", HomotopyEquivalence::is_connected(&complex));
println!("β₀ (components): {}", complex.betti_0());
```

### Radius Sweep (Filtration)

```rust
// Build complexes at multiple radii to see topology evolve
let sweep = RadiusSweep::sweep(&points, &[0.1, 0.5, 1.0, 2.0], 2);

for (i, complex) in sweep.iter().enumerate() {
    println!("r={:.1}: {} simplices, {} components",
        [0.1, 0.5, 1.0, 2.0][i],
        complex.num_simplices(),
        complex.betti_0(),
    );
}

// Find the connection radius (minimum r for a connected complex)
let r_conn = RadiusSweep::connection_radius(&points, 2);
println!("Connection radius: {:.3}", r_conn);
```

### Ball Intersection Checks

```rust
use cech_complex::BallIntersection;

// Check if two specific balls intersect at radius r
let points = vec![vec![0.0], vec![1.5], vec![5.0]];
let r = 1.0;

println!("Balls 0,1 intersect: {}", BallIntersection::pair_intersects(&points, 0, 1, r));
println!("Balls 1,2 intersect: {}", BallIntersection::pair_intersects(&points, 1, 2, r));

// Check triple intersection
println!("Triple intersection: {}",
    BallIntersection::intersects(&points, &[0, 1, 2], r));
```

## API Reference

### BallIntersection

| Method | Returns | Description |
|--------|---------|-------------|
| `BallIntersection::pair_intersects(points, i, j, r)` | `bool` | Do balls around points i, j intersect at radius r? |
| `BallIntersection::intersects(points, indices, r)` | `bool` | Do all balls around given indices intersect? |

### Nerve

| Method | Returns | Description |
|--------|---------|-------------|
| `Nerve::build(n, max_dim, predicate)` | `HashSet<Vec<usize>>` | Build nerve from an arbitrary intersection predicate |

### CechComplex

| Method | Returns | Description |
|--------|---------|-------------|
| `CechComplex::build(points, r, max_dim)` | `CechComplex` | Build the Čech complex at radius r |
| `complex.num_vertices()` | `usize` | Number of input points |
| `complex.num_simplices()` | `usize` | Total number of simplices |
| `complex.simplices_of_dim(d)` | `Vec<&Vec<usize>>` | All d-dimensional simplices |
| `complex.has_simplex(vertices)` | `bool` | Check if simplex exists |
| `complex.betti_0()` | `usize` | Number of connected components |
| `complex.radius` | `f64` | The radius used to build this complex |

### HomotopyEquivalence

| Method | Returns | Description |
|--------|---------|-------------|
| `HomotopyEquivalence::nerve_theorem_holds(dim)` | `bool` | Always true for balls in ℝᵈ (convex) |
| `HomotopyEquivalence::is_connected(complex)` | `bool` | Is the complex connected? |

### RadiusSweep

| Method | Returns | Description |
|--------|---------|-------------|
| `RadiusSweep::sweep(points, radii, max_dim)` | `Vec<CechComplex>` | Build complexes at multiple radii |
| `RadiusSweep::connection_radius(points, max_dim)` | `f64` | Minimum radius for connectivity |
| `RadiusSweep::simplex_counts(sweep)` | `Vec<usize>` | Simplex counts at each radius |

## Mathematical Background

### The Čech Complex

Given a finite point set P = {p₁, ..., pₙ} ⊂ ℝᵈ and radius r ≥ 0, the Čech complex Čech(r) has:
- A **vertex** for each point pᵢ
- A **k-simplex** {pᵢ₀, ..., pᵢₖ} when ∩ⱼ₌₀ᵏ B(pᵢⱼ, r) ≠ ∅

where B(p, r) = {x ∈ ℝᵈ : ||x − p|| ≤ r} is the closed ball of radius r centered at p.

### The Nerve Theorem

**Theorem** (Borsuk, 1948; refined by Hatcher): If 𝒰 is a finite collection of convex sets in ℝᵈ, then the nerve N(𝒰) is homotopy equivalent to the union ∪𝒰.

Since balls in Euclidean space are convex, and intersections of convex sets are convex, the Čech complex is homotopy equivalent to the union of balls:
```
Čech(r) ≃ ∪ᵢ B(pᵢ, r)
```

This means the Čech complex captures the **exact topology** (homotopy type) of the ball union — no approximation error.

### Čech vs. Vietoris-Rips

| Property | Čech | Vietoris-Rips |
|----------|------|---------------|
| Condition | All-way ball intersection | Pairwise ball intersection |
| Topology | Exact (Nerve Theorem) | Approximation |
| Size | Smaller | Larger (more simplices) |
| Computation | Harder (circumradius check) | Easier (pairwise distance) |
| Relationship | Čech(r) ⊂ VR(r) ⊂ Čech(r·√(2d/k)) | |

### Circumradius Criterion

For k+1 points to form a k-simplex in the Čech complex, their balls must all intersect. This is equivalent to the circumradius of the points being ≤ r. For 2 points (an edge), this is simply: distance ≤ 2r. For 3+ points, it requires computing the circumradius of the simplex.

## Installation

```bash
cargo add cech-complex
```

Or add to your `Cargo.toml`:

```toml
[dependencies]
cech-complex = "0.1.0"
```

## Related Crates

- [`witness-complex`](https://github.com/SuperInstance/witness-complex) — Scalable approximation via landmarks
- [`mapper-graph`](https://github.com/SuperInstance/mapper-graph) — Mapper algorithm for topological summaries
- [`persistence-landscape`](https://github.com/SuperInstance/persistence-landscape) — Persistence landscapes for statistical TDA
- [`betti-curve`](https://github.com/SuperInstance/betti-curve) — Betti curves and Euler characteristic curves

## License

MIT © [SuperInstance](https://github.com/SuperInstance)

---

*Part of the [Exocortex](https://github.com/SuperInstance/exocortex) project — persistent cognitive substrate for multi-agent systems.*
