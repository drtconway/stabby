use std::{
    collections::{HashMap, VecDeque},
    fmt::Display,
};

use ransel::{rank::Rank, select::Select, sparse::Sparse};

use crate::dense::{DenseInterval, DenseStabby};

/// The Interval struct represents a closed interval on an (unsigned) integer domain.
#[derive(Clone, Copy, Eq, PartialOrd, Ord, Default, Hash, PartialEq, Debug)]
pub struct Interval {
    pub first: u64,
    pub last: u64,
}

impl Interval {
    pub fn new(first: u64, last: u64) -> Interval {
        debug_assert!(first <= last);
        Interval { first, last }
    }

    pub fn zero() -> Interval {
        Interval { first: 0, last: 0 }
    }
}

impl Display for Interval {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "[{}, {}]", self.first, self.last)
    }
}

pub struct Stabby {
    domain: Sparse,
    dense: DenseStabby,
}

///
impl Stabby {
    pub fn new(xs: &[Interval]) -> Stabby {
        let domain = Self::make_domain(xs);
        let mut ys: Vec<DenseInterval> = Vec::new();
        let mut y_max = 0;
        for x in xs.iter() {
            let (y_f, y_l) = domain.rank_2(x.first as u64, x.last as u64);
            let y = DenseInterval::new(y_f * 2, y_l * 2);
            if y.last > y_max {
                y_max = y.last;
            }
            ys.push(y);
        }
        ys.sort();
        let dense = DenseStabby::new(y_max as usize + 1, &ys);

        Stabby { domain, dense }
    }

    fn make_domain(xs: &[Interval]) -> Sparse {
        let mut ys: Vec<u64> = Vec::new();
        ys.push(0);
        for x in xs.iter() {
            ys.push(x.first);
            ys.push(x.last);
        }
        ys.sort();
        ys.dedup();
        let x = 1 + *ys.last().unwrap();
        let b = 1 + (64 - x.leading_zeros()) as usize;
        let n = ys.len();
        Sparse::new(b, n, &ys)
    }

    pub fn stabs(&self, q: u64) -> bool {
        let qd = self.sparse_to_dense(q);
        self.dense.stabs(qd)
    }

    pub fn stab(&self, q: u64) -> Vec<Interval> {
        let qd = self.sparse_to_dense(q);
        let ys = self.dense.stab(qd);
        let mut xs: Vec<Interval> = Vec::new();
        for y in ys.iter() {
            let y1 = self.dense_to_sparse(y.first);
            let y2 = self.dense_to_sparse(y.last);
            xs.push(Interval::new(y1, y2));
        }
        xs
    }

    pub fn stab_interval(&self, q: &Interval) -> Vec<Interval> {
        let qd = DenseInterval::new(self.sparse_to_dense(q.first), self.sparse_to_dense(q.last));
        println!("qd = {:?}", qd);
        println!("dense = {:?}", self.dense);
        let ys = self.dense.stab_interval(&qd);
        let mut xs: Vec<Interval> = Vec::new();
        for y in ys.iter() {
            let y1 = self.dense_to_sparse(y.first as usize);
            let y2 = self.dense_to_sparse(y.last as usize);
            xs.push(Interval::new(y1, y2));
        }
        xs
    }

    fn sparse_to_dense(&self, x: u64) -> usize {
        let (r1, r2) = self.domain.rank_2(x as u64, (x + 1) as u64);
        let q = if r1 != r2 { r1 * 2 } else { r1 * 2 - 1 };
        q
    }

    fn dense_to_sparse(&self, x: usize) -> u64 {
        self.domain.select(x / 2)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_stabby_1() {
        let mutyh = vec![
            Interval::new(45329163, 45329437),
            Interval::new(45330516, 45330557),
            Interval::new(45331182, 45331334),
            Interval::new(45331420, 45331556),
            Interval::new(45332763, 45332834),
            Interval::new(45334391, 45334511),
            Interval::new(45340219, 45340447),
        ];
        let idx = Stabby::new(&mutyh);
        assert!(idx.stabs(45_331_258));
        assert!(!idx.stabs(45_332_840));
        assert_eq!(
            idx.stab(45_331_258),
            vec![Interval::new(45_331_182, 45_331_334)]
        );
        assert_eq!(
            idx.stab_interval(&Interval::new(45_331_151, 45_331_880)),
            vec![
                Interval::new(45_331_182, 45_331_334),
                Interval::new(45_331_420, 45_331_556)
            ]
        );
    }

    #[test]
    fn test_stabby_4() {
        let src: Vec<Interval> = Vec::from([
            Interval::new(55519604, 55519916),
            Interval::new(55519604, 55545079),
            Interval::new(55519646, 55519916),
            Interval::new(55519646, 55539543),
            Interval::new(55519690, 55519916),
            Interval::new(55519690, 55545079),
            Interval::new(55519718, 55519916),
            Interval::new(55519718, 55545079),
            Interval::new(55519731, 55519916),
            Interval::new(55519731, 55535914),
            Interval::new(55520405, 55520479),
            Interval::new(55522102, 55522166),
            Interval::new(55523721, 55523822),
            Interval::new(55528878, 55528992),
            Interval::new(55533873, 55533960),
            Interval::new(55535712, 55535763),
            Interval::new(55535712, 55535914),
            Interval::new(55537483, 55537585),
            Interval::new(55538765, 55539543),
            Interval::new(55543938, 55544074),
            Interval::new(55544025, 55544074),
            Interval::new(55544220, 55544369),
            Interval::new(55544907, 55545079),
            Interval::new(55547292, 55550006),
            Interval::new(55547292, 55617614),
            Interval::new(55547292, 55617622),
            Interval::new(55547292, 55618880),
            Interval::new(55548378, 55550006),
            Interval::new(55548378, 55617608),
            Interval::new(55558775, 55558968),
            Interval::new(55564313, 55564497),
            Interval::new(55564902, 55565041),
            Interval::new(55565703, 55565850),
            Interval::new(55568194, 55568363),
            Interval::new(55573619, 55573777),
            Interval::new(55577315, 55577356),
            Interval::new(55578247, 55578342),
            Interval::new(55579679, 55579781),
            Interval::new(55581567, 55581698),
            Interval::new(55585051, 55585167),
            Interval::new(55586618, 55586734),
            Interval::new(55588879, 55588956),
            Interval::new(55598416, 55599039),
            Interval::new(55603978, 55604076),
            Interval::new(55615451, 55615506),
            Interval::new(55617144, 55617608),
            Interval::new(55617144, 55617614),
            Interval::new(55617144, 55617622),
            Interval::new(55617869, 55618444),
            Interval::new(55634061, 55636392),
            Interval::new(55634061, 55693844),
            Interval::new(55634061, 55693863),
            Interval::new(55637552, 55637599),
            Interval::new(55640627, 55640705),
            Interval::new(55643158, 55643213),
            Interval::new(55643319, 55643425),
            Interval::new(55644637, 55644720),
            Interval::new(55645349, 55645432),
            Interval::new(55646259, 55646322),
            Interval::new(55646415, 55646486),
            Interval::new(55647347, 55647453),
            Interval::new(55654900, 55654953),
            Interval::new(55656131, 55656220),
            Interval::new(55656305, 55656371),
            Interval::new(55660157, 55660193),
            Interval::new(55661956, 55662026),
            Interval::new(55666991, 55667093),
            Interval::new(55667862, 55667958),
            Interval::new(55671319, 55671376),
            Interval::new(55671995, 55672046),
            Interval::new(55672893, 55673079),
            Interval::new(55679682, 55679795),
            Interval::new(55680712, 55680759),
            Interval::new(55680855, 55680918),
            Interval::new(55683785, 55683834),
            Interval::new(55684943, 55685048),
            Interval::new(55684994, 55685048),
            Interval::new(55684994, 55693823),
            Interval::new(55686370, 55686444),
            Interval::new(55687645, 55687705),
            Interval::new(55693663, 55693823),
            Interval::new(55693663, 55693844),
        ]);
        let s: Stabby = Stabby::new(&src);
        let mut q = 55519000;
        while q < 55700000 {
            let mut expected: Vec<Interval> = Vec::new();
            for ivl in src.iter() {
                if ivl.first <= q && q <= ivl.last {
                    expected.push(*ivl);
                }
            }
            assert_eq!(s.stab(q), expected);
            q += 150;
        }
    }
}
