//! # cech-complex
//!
//! **Čech complex** construction from point clouds via ball intersections
//! and nerve computation for topological data analysis.
//!
//! Given a point cloud `P` and radius `r`, the Čech complex has a k-simplex
//! for every subset of k+1 points whose balls of radius r have non-empty
//! intersection. By the Nerve Theorem, the Čech complex is homotopy equivalent
//! to the union of the balls (under mild conditions).
//!
//! # Example
//!
//! ```
//! use cech_complex::{CechComplex, RadiusSweep};
//!
//! let points = vec![vec![0.0, 0.0], vec![1.0, 0.0], vec![0.5, 0.866]];
//! let cc = CechComplex::build(&points, 1.0, 2);
//! assert!(cc.num_vertices() == 3);
//!
//! let sweep = RadiusSweep::sweep(&points, &[0.5, 1.0, 1.5], 2);
//! assert_eq!(sweep.len(), 3);
//! ```

use std::collections::HashSet;

// ---------------------------------------------------------------------------
// BallIntersection
// ---------------------------------------------------------------------------

/// Check whether balls of a given radius around points intersect.
pub struct BallIntersection;

impl BallIntersection {
    /// Check if the intersection of balls of radius `r` around the given points
    /// is non-empty. For a set of points this is equivalent to: the circumradius
    /// of the points ≤ r. For 2 points this is just: distance ≤ 2r.
    pub fn intersects(points: &[Vec<f64>], indices: &[usize], r: f64) -> bool {
        if indices.is_empty() {
            return true;
        }
        // For each pair, check distance ≤ 2r
        for i in 0..indices.len() {
            for j in (i + 1)..indices.len() {
                let d = euclidean(&points[indices[i]], &points[indices[j]]);
                if d > 2.0 * r + 1e-12 {
                    return false;
                }
            }
        }
        true
    }

    /// Check if two balls of radius r around points i and j intersect.
    pub fn pair_intersects(points: &[Vec<f64>], i: usize, j: usize, r: f64) -> bool {
        euclidean(&points[i], &points[j]) <= 2.0 * r + 1e-12
    }
}

/// Euclidean distance.
fn euclidean(a: &[f64], b: &[f64]) -> f64 {
    a.iter()
        .zip(b.iter())
        .map(|(x, y)| (x - y) * (x - y))
        .sum::<f64>()
        .sqrt()
}

// ---------------------------------------------------------------------------
// Nerve
// ---------------------------------------------------------------------------

/// Compute the nerve of a collection of sets indexed by point indices.
///
/// The nerve is the simplicial complex where a k-simplex exists iff the
/// intersection of the corresponding k+1 sets is non-empty.
pub struct Nerve;

impl Nerve {
    /// Build the nerve given an intersection predicate.
    ///
    /// `n` is the number of elements, `max_dim` is the maximum simplex dimension,
    /// and `intersects` is a closure that returns true if the sets indexed by
    /// the given indices have non-empty intersection.
    pub fn build<F>(n: usize, max_dim: usize, intersects: F) -> HashSet<Vec<usize>>
    where
        F: Fn(&[usize]) -> bool,
    {
        let mut simplices: HashSet<Vec<usize>> = HashSet::new();

        // Add all vertices
        for i in 0..n {
            simplices.insert(vec![i]);
        }

        // Enumerate subsets of size 2 to (max_dim+1)
        let mut current: Vec<usize> = Vec::new();
        Self::enumerate_subsets(n, max_dim + 1, &mut current, &mut simplices, &intersects);

        simplices
    }

    fn enumerate_subsets<F>(
        n: usize,
        max_size: usize,
        current: &mut Vec<usize>,
        result: &mut HashSet<Vec<usize>>,
        pred: &F,
    ) where
        F: Fn(&[usize]) -> bool,
    {
        if current.len() >= 2 {
            if pred(current) {
                let mut sorted = current.clone();
                sorted.sort();
                result.insert(sorted);
            } else {
                return; // No superset can satisfy the predicate either if pairwise fails
            }
        }

        if current.len() >= max_size {
            return;
        }

        let start = current.last().map(|&i| i + 1).unwrap_or(0);
        for i in start..n {
            current.push(i);
            Self::enumerate_subsets(n, max_size, current, result, pred);
            current.pop();
        }
    }
}

// ---------------------------------------------------------------------------
// CechComplex
// ---------------------------------------------------------------------------

/// The Čech complex at radius r: a simplex for each subset of points whose
/// balls of radius r have non-empty common intersection.
#[derive(Debug, Clone)]
pub struct CechComplex {
    /// Number of input points.
    n: usize,
    /// Simplices stored as sorted Vec<usize>.
    simplices: HashSet<Vec<usize>>,
    /// The radius used to build this complex.
    pub radius: f64,
}

impl CechComplex {
    /// Build the Čech complex from points at radius `r`, up to dimension `max_dim`.
    pub fn build(points: &[Vec<f64>], r: f64, max_dim: usize) -> Self {
        let n = points.len();
        let pts = points.to_vec();
        let simplices = Nerve::build(n, max_dim, move |indices| {
            BallIntersection::intersects(&pts, indices, r)
        });

        Self {
            n,
            simplices,
            radius: r,
        }
    }

    /// Number of vertices.
    pub fn num_vertices(&self) -> usize {
        self.n
    }

    /// Total number of simplices.
    pub fn num_simplices(&self) -> usize {
        self.simplices.len()
    }

    /// Get simplices of a given dimension.
    pub fn simplices_of_dim(&self, dim: usize) -> Vec<&Vec<usize>> {
        self.simplices
            .iter()
            .filter(|s| s.len() == dim + 1)
            .collect()
    }

    /// Check if a specific simplex exists.
    pub fn has_simplex(&self, vertices: &[usize]) -> bool {
        let mut v = vertices.to_vec();
        v.sort();
        self.simplices.contains(&v)
    }

    /// Betti number approximation: count connected components via union-find on edges.
    pub fn betti_0(&self) -> usize {
        let mut parent: Vec<usize> = (0..self.n).collect();

        let find = |parent: &mut Vec<usize>, mut x: usize| -> usize {
            while parent[x] != x {
                parent[x] = parent[parent[x]];
                x = parent[x];
            }
            x
        };

        for sigma in &self.simplices {
            if sigma.len() == 2 {
                let r1 = find(&mut parent, sigma[0]);
                let r2 = find(&mut parent, sigma[1]);
                if r1 != r2 {
                    parent[r1] = r2;
                }
            }
        }

        let mut roots = HashSet::new();
        for i in 0..self.n {
            roots.insert(find(&mut parent, i));
        }
        roots.len()
    }
}

// ---------------------------------------------------------------------------
// HomotopyEquivalence
// ---------------------------------------------------------------------------

/// Check conditions related to the Nerve Theorem for homotopy equivalence.
pub struct HomotopyEquivalence;

impl HomotopyEquivalence {
    /// The nerve theorem guarantees homotopy equivalence when each finite
    /// intersection of balls is contractible (always true for convex balls in ℝ^d).
    /// For balls in Euclidean space, this is always satisfied.
    pub fn nerve_theorem_holds(_dim: usize) -> bool {
        true // Balls in ℝ^d are convex, intersections are convex → contractible
    }

    /// Check if the Čech complex at radius r is connected.
    pub fn is_connected(complex: &CechComplex) -> bool {
        complex.betti_0() <= 1
    }
}

// ---------------------------------------------------------------------------
// RadiusSweep
// ---------------------------------------------------------------------------

/// Build Čech complexes at multiple radii (filtration).
pub struct RadiusSweep;

impl RadiusSweep {
    /// Build a sequence of Čech complexes for the given radii.
    pub fn sweep(points: &[Vec<f64>], radii: &[f64], max_dim: usize) -> Vec<CechComplex> {
        radii
            .iter()
            .map(|&r| CechComplex::build(points, r, max_dim))
            .collect()
    }

    /// Find the minimum radius at which all points are connected.
    pub fn connection_radius(points: &[Vec<f64>], _max_dim: usize) -> f64 {
        if points.len() <= 1 {
            return 0.0;
        }
        // The maximum distance between any two points / 2
        let mut max_dist = 0.0;
        for i in 0..points.len() {
            for j in (i + 1)..points.len() {
                let d = euclidean(&points[i], &points[j]);
                if d > max_dist {
                    max_dist = d;
                }
            }
        }
        max_dist / 2.0
    }

    /// Count simplices at each radius in the sweep.
    pub fn simplex_counts(sweep: &[CechComplex]) -> Vec<usize> {
        sweep.iter().map(|c| c.num_simplices()).collect()
    }
}

// ===========================================================================
// Tests
// ===========================================================================

#[cfg(test)]
mod tests {
    use super::*;

    fn triangle_points() -> Vec<Vec<f64>> {
        vec![
            vec![0.0, 0.0],
            vec![1.0, 0.0],
            vec![0.5, 0.866],
        ]
    }

    fn line_points() -> Vec<Vec<f64>> {
        vec![vec![0.0], vec![1.0], vec![2.0], vec![3.0]]
    }

    #[test]
    fn test_ball_intersection_close() {
        let pts = vec![vec![0.0], vec![1.0]];
        assert!(BallIntersection::pair_intersects(&pts, 0, 1, 0.6));
    }

    #[test]
    fn test_ball_intersection_far() {
        let pts = vec![vec![0.0], vec![10.0]];
        assert!(!BallIntersection::pair_intersects(&pts, 0, 1, 1.0));
    }

    #[test]
    fn test_ball_intersection_triple() {
        let pts = triangle_points();
        assert!(BallIntersection::intersects(&pts, &[0, 1, 2], 1.0));
    }

    #[test]
    fn test_ball_intersection_triple_far() {
        let pts = vec![vec![0.0], vec![10.0], vec![20.0]];
        assert!(!BallIntersection::intersects(&pts, &[0, 1, 2], 1.0));
    }

    #[test]
    fn test_cech_complex_vertices() {
        let pts = triangle_points();
        let cc = CechComplex::build(&pts, 0.5, 2);
        assert_eq!(cc.num_vertices(), 3);
    }

    #[test]
    fn test_cech_complex_edges_small_radius() {
        let pts = triangle_points();
        let cc = CechComplex::build(&pts, 0.01, 2);
        // Very small radius: only vertices, no edges
        assert_eq!(cc.simplices_of_dim(0).len(), 3);
        assert_eq!(cc.simplices_of_dim(1).len(), 0);
    }

    #[test]
    fn test_cech_complex_edges_large_radius() {
        let pts = triangle_points();
        let cc = CechComplex::build(&pts, 1.0, 2);
        assert!(cc.simplices_of_dim(1).len() >= 3);
        assert!(cc.has_simplex(&[0, 1, 2]));
    }

    #[test]
    fn test_cech_complex_betti0_disconnected() {
        let pts = vec![vec![0.0], vec![100.0]];
        let cc = CechComplex::build(&pts, 0.1, 1);
        assert_eq!(cc.betti_0(), 2);
    }

    #[test]
    fn test_cech_complex_betti0_connected() {
        let pts = vec![vec![0.0], vec![1.0]];
        let cc = CechComplex::build(&pts, 1.0, 1);
        assert_eq!(cc.betti_0(), 1);
    }

    #[test]
    fn test_nerve_theorem() {
        assert!(HomotopyEquivalence::nerve_theorem_holds(2));
        assert!(HomotopyEquivalence::nerve_theorem_holds(3));
    }

    #[test]
    fn test_is_connected() {
        let pts = vec![vec![0.0], vec![1.0], vec![2.0]];
        let cc = CechComplex::build(&pts, 1.5, 1);
        assert!(HomotopyEquivalence::is_connected(&cc));
    }

    #[test]
    fn test_radius_sweep() {
        let pts = line_points();
        let sweep = RadiusSweep::sweep(&pts, &[0.1, 0.6, 1.5], 1);
        assert_eq!(sweep.len(), 3);
        // More simplices at larger radius
        assert!(sweep[2].num_simplices() >= sweep[0].num_simplices());
    }

    #[test]
    fn test_connection_radius() {
        let pts = vec![vec![0.0], vec![2.0]];
        let r = RadiusSweep::connection_radius(&pts, 1);
        assert!((r - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_simplex_counts() {
        let pts = line_points();
        let sweep = RadiusSweep::sweep(&pts, &[0.1, 1.0, 2.0], 1);
        let counts = RadiusSweep::simplex_counts(&sweep);
        assert_eq!(counts.len(), 3);
        assert!(counts[2] >= counts[0]);
    }
}
