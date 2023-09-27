use std::{
    collections::{HashMap, VecDeque},
    hash::Hash,
};

use crate::listy::{Listy, ListyElement};

#[derive(Clone, Copy, Eq, PartialOrd, Ord, Default, Hash, PartialEq, Debug)]
pub struct DenseInterval {
    pub first: usize,
    pub last: usize,
}

impl DenseInterval {
    pub fn new(first: usize, last: usize) -> DenseInterval {
        DenseInterval { first, last }
    }

    pub fn zero() -> DenseInterval {
        DenseInterval { first: 0, last: 0 }
    }
}

/// Internal to the Schmidt algorithm, we compute the "smaller" relationship.
///
fn make_smaller(
    items: &[DenseInterval],
) -> (
    Vec<DenseInterval>,
    HashMap<DenseInterval, Vec<DenseInterval>>,
) {
    let n = items.len();
    let mut i: usize = 0;
    let mut basic: Vec<DenseInterval> = Vec::new();
    let mut smaller: HashMap<DenseInterval, Vec<DenseInterval>> = HashMap::new();
    while i < n {
        let mut j: usize = i;
        while j + 1 < n && items[j].first == items[j + 1].first {
            j += 1;
        }
        if j > i {
            let v = &items[i..j];
            if v.len() > 0 {
                smaller.insert(items[j], Vec::from(v));
            }
        }
        basic.push(items[j]);
        i = j + 1;
    }
    (basic, smaller)
}

/// The DenseStabby data structure represents the set of intervals over the dense domain.
///
#[derive(Debug)]
pub struct DenseStabby {
    smaller: HashMap<DenseInterval, Vec<DenseInterval>>,
    start: HashMap<usize, DenseInterval>,
    start2: HashMap<usize, DenseInterval>,
    parent: HashMap<DenseInterval, DenseInterval>,
    last: HashMap<DenseInterval, DenseInterval>,
    left: HashMap<DenseInterval, DenseInterval>,
}

impl DenseStabby {
    pub fn new(q_max: usize, items: &[DenseInterval]) -> DenseStabby {
        let (basic, smaller) = make_smaller(items);
        let mut event: Vec<Vec<DenseInterval>> = Vec::new();
        event.resize_with(q_max + 1, || Vec::new());

        for item in basic.iter() {
            println!("queuing {}, {} ({})", item.first, item.last, q_max);
            event[item.last].push(*item);
            event[item.first].push(*item);
        }

        let mut start: HashMap<usize, DenseInterval> = HashMap::new();
        let mut startx: Vec<Option<DenseInterval>> = Vec::new();
        startx.resize(q_max + 1, None);

        let mut start2: HashMap<usize, DenseInterval> = HashMap::new();
        let mut start2x: Vec<Option<DenseInterval>> = Vec::new();
        start2x.resize(q_max + 1, None);

        let mut parent: HashMap<DenseInterval, DenseInterval> = HashMap::new();
        let mut last: HashMap<DenseInterval, DenseInterval> = HashMap::new();
        let mut left: HashMap<DenseInterval, DenseInterval> = HashMap::new();

        let mut l: Listy<DenseInterval> = Listy::new();
        let mut saved: HashMap<DenseInterval, ListyElement<DenseInterval>> = HashMap::new();
        let mut rml: usize = 0;

        for q in 0..=q_max {
            match l.back() {
                None => {}
                Some(a) => {
                    start.insert(q, *a);
                    startx[q] = Some(*a);
                }
            }
            while let Some(a) = event[q].pop() {
                if a.first == q {
                    start.insert(q, a);
                    startx[q] = Some(a);
                    let ptr = l.push_back(a);
                    saved.insert(a, ptr);
                } else {
                    let mut p = DenseInterval::zero();
                    if let Some(a_ptr) = saved.get(&a) {
                        if let Some(b_ptr) = l.prev(a_ptr) {
                            let b = *l.get(&b_ptr);
                            p = b;
                        }
                    }

                    parent.insert(a, p);
                    if let Some(q) = last.get(&p) {
                        left.insert(a, *q);
                    }
                    last.insert(p, a);

                    if let Some(a_ptr) = saved.get(&a) {
                        l.remove(&a_ptr);
                    }

                    saved.remove(&a);
                }
            }
            while rml + 1 < basic.len() && basic[rml + 1].first <= q {
                rml += 1;
            }
            if basic[rml].first <= q {
                start2.insert(q, basic[rml]);
                start2x[q] = Some(basic[rml]);
            }
        }

        println!("start = {:?}", startx);
        println!("start2 = {:?}", start2x);

        DenseStabby {
            smaller: smaller,
            start: start,
            start2: start2,
            parent: parent,
            last: last,
            left: left,
        }
    }

    /// A quick test to see if a position is included in any intervals
    /// without determining which specific intervals.
    pub fn stabs(&self, q: usize) -> bool {
        self.start.contains_key(&q)
    }

    pub fn stab(&self, q: usize) -> Vec<DenseInterval> {
        let mut res: Vec<DenseInterval> = Vec::new();
        let mut kew: VecDeque<DenseInterval> = VecDeque::new();
        let mut ov: Option<&DenseInterval> = self.start.get(&q);
        loop {
            match ov {
                None => {
                    break;
                }
                Some(v) => {
                    if *v == DenseInterval::zero() {
                        break;
                    }
                    kew.push_front(*v);
                    ov = self.parent.get(v);
                }
            }
        }
        loop {
            match kew.pop_back() {
                None => {
                    break;
                }
                Some(a) => {
                    res.push(a);
                    match self.smaller.get(&a) {
                        None => {}
                        Some(s) => {
                            for r in s.iter().rev() {
                                if r.last < q {
                                    break;
                                }
                                res.push(*r);
                            }
                        }
                    }
                    let mut ot = self.left.get(&a);
                    loop {
                        match ot {
                            None => {
                                break;
                            }
                            Some(t) => {
                                if t.last < q {
                                    break;
                                }
                                kew.push_back(*t);
                                ot = self.last.get(t);
                            }
                        }
                    }
                }
            }
        }
        res.reverse();
        res
    }

    pub fn stab_interval(&self, qi: &DenseInterval) -> Vec<DenseInterval> {
        let lq = qi.first;
        let rq = qi.last;

        let mut ot: Option<&DenseInterval> = None;
        match self.start.get(&lq) {
            None => {}
            Some(u) => {
                ot = Some(u);
            }
        }
        match self.start2.get(&rq) {
            None => {}
            Some(u) => match ot {
                None => {
                    ot = Some(u);
                }
                Some(t) => {
                    if t.first < u.first {
                        ot = Some(u);
                    }
                }
            },
        }

        let mut res: Vec<DenseInterval> = Vec::new();
        match ot {
            None => {
                return res;
            }
            Some(t) => {
                if t.last < lq {
                    return res;
                }
            }
        }

        let mut kew: VecDeque<DenseInterval> = VecDeque::new();
        loop {
            match ot {
                None => {
                    break;
                }
                Some(t) => {
                    if *t == DenseInterval::zero() {
                        break;
                    }
                    kew.push_front(*t);
                    ot = self.parent.get(&t);
                }
            }
        }

        loop {
            match kew.pop_back() {
                None => {
                    break;
                }
                Some(a) => {
                    res.push(a);

                    match self.smaller.get(&a) {
                        None => {}
                        Some(s) => {
                            for r in s.iter().rev() {
                                if r.last < lq {
                                    break;
                                }
                                res.push(*r);
                            }
                        }
                    }

                    ot = self.left.get(&a);

                    loop {
                        match ot {
                            None => {
                                break;
                            }
                            Some(t) => {
                                if t.last < lq {
                                    break;
                                }
                                kew.push_back(*t);
                                ot = self.last.get(t);
                            }
                        }
                    }
                }
            }
        }
        res.reverse();
        res
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_make_smaller() {
        let items: [DenseInterval; 3] = [
            DenseInterval { first: 1, last: 3 },
            DenseInterval { first: 2, last: 5 },
            DenseInterval { first: 2, last: 7 },
        ];
        let (basic, _smaller) = make_smaller(&items);
        assert_eq!(basic.len(), 2);
    }

    #[test]
    fn test_stabby_1() {
        let src: Vec<DenseInterval> = Vec::from([
            DenseInterval::new(1, 2),
            DenseInterval::new(1, 4),
            DenseInterval::new(3, 5),
        ]);
        let s = DenseStabby::new(6, &src);
        let w = DenseInterval::new(1, 4);
        let v = Vec::from([DenseInterval::new(1, 2)]);
        assert_eq!(s.smaller.len(), 1);
        assert_eq!(s.smaller.get(&w), Some(&v));
        assert_eq!(s.stabs(3), true);
        assert_eq!(s.stabs(0), false);
        assert_eq!(
            s.stab(3),
            Vec::from([DenseInterval::new(1, 4), DenseInterval::new(3, 5)])
        );
    }

    #[test]
    fn test_stabby_2() {
        let mut src: Vec<DenseInterval> = Vec::from([
            DenseInterval::new(7, 62),
            DenseInterval::new(38, 57),
            DenseInterval::new(22, 65),
            DenseInterval::new(5, 88),
            DenseInterval::new(23, 43),
            DenseInterval::new(73, 75),
            DenseInterval::new(36, 89),
            DenseInterval::new(31, 55),
            DenseInterval::new(17, 91),
            DenseInterval::new(24, 80),
            DenseInterval::new(40, 47),
            DenseInterval::new(52, 67),
            DenseInterval::new(16, 45),
            DenseInterval::new(86, 90),
            DenseInterval::new(46, 56),
            DenseInterval::new(33, 84),
            DenseInterval::new(21, 56),
            DenseInterval::new(15, 93),
            DenseInterval::new(48, 51),
            DenseInterval::new(68, 96),
            DenseInterval::new(35, 64),
            DenseInterval::new(63, 78),
            DenseInterval::new(71, 99),
            DenseInterval::new(3, 4),
            DenseInterval::new(18, 85),
            DenseInterval::new(37, 83),
            DenseInterval::new(10, 25),
            DenseInterval::new(32, 66),
            DenseInterval::new(9, 49),
            DenseInterval::new(3, 13),
            DenseInterval::new(6, 94),
            DenseInterval::new(75, 98),
            DenseInterval::new(8, 69),
            DenseInterval::new(39, 74),
            DenseInterval::new(11, 41),
            DenseInterval::new(14, 59),
            DenseInterval::new(24, 75),
            DenseInterval::new(1, 2),
            DenseInterval::new(50, 53),
            DenseInterval::new(0, 42),
            DenseInterval::new(49, 76),
            DenseInterval::new(25, 55),
            DenseInterval::new(60, 69),
            DenseInterval::new(15, 18),
            DenseInterval::new(52, 79),
            DenseInterval::new(15, 74),
            DenseInterval::new(30, 82),
            DenseInterval::new(34, 81),
            DenseInterval::new(39, 97),
            DenseInterval::new(51, 81),
            DenseInterval::new(29, 87),
            DenseInterval::new(61, 92),
            DenseInterval::new(12, 70),
            DenseInterval::new(26, 54),
            DenseInterval::new(12, 52),
            DenseInterval::new(52, 54),
            DenseInterval::new(12, 72),
            DenseInterval::new(12, 58),
            DenseInterval::new(44, 88),
        ]);
        src.sort();
        let s = DenseStabby::new(100, &src);
        for q in 0..=100 {
            let mut expected: Vec<DenseInterval> = Vec::new();
            for ivl in src.iter() {
                if ivl.first <= q && q <= ivl.last {
                    expected.push(*ivl);
                }
            }
            assert_eq!(s.stabs(q), expected.len() > 0);
            assert_eq!(s.stab(q), expected);
        }
    }

    #[test]
    fn test_stabby_3() {
        let mut src: Vec<DenseInterval> = Vec::from([
            DenseInterval::new(7, 62),
            DenseInterval::new(38, 57),
            DenseInterval::new(22, 65),
            DenseInterval::new(5, 88),
            DenseInterval::new(23, 43),
            DenseInterval::new(73, 75),
            DenseInterval::new(36, 89),
            DenseInterval::new(31, 55),
            DenseInterval::new(17, 91),
            DenseInterval::new(24, 80),
            DenseInterval::new(40, 47),
            DenseInterval::new(52, 67),
            DenseInterval::new(16, 45),
            DenseInterval::new(86, 90),
            DenseInterval::new(46, 56),
            DenseInterval::new(33, 84),
            DenseInterval::new(21, 56),
            DenseInterval::new(15, 93),
            DenseInterval::new(48, 51),
            DenseInterval::new(68, 96),
            DenseInterval::new(35, 64),
            DenseInterval::new(63, 78),
            DenseInterval::new(71, 99),
            DenseInterval::new(3, 4),
            DenseInterval::new(18, 85),
            DenseInterval::new(37, 83),
            DenseInterval::new(10, 25),
            DenseInterval::new(32, 66),
            DenseInterval::new(9, 49),
            DenseInterval::new(3, 13),
            DenseInterval::new(6, 94),
            DenseInterval::new(75, 98),
            DenseInterval::new(8, 69),
            DenseInterval::new(39, 74),
            DenseInterval::new(11, 41),
            DenseInterval::new(14, 59),
            DenseInterval::new(24, 75),
            DenseInterval::new(1, 2),
            DenseInterval::new(50, 53),
            DenseInterval::new(0, 42),
            DenseInterval::new(49, 76),
            DenseInterval::new(25, 55),
            DenseInterval::new(60, 69),
            DenseInterval::new(15, 18),
            DenseInterval::new(52, 79),
            DenseInterval::new(15, 74),
            DenseInterval::new(30, 82),
            DenseInterval::new(34, 81),
            DenseInterval::new(39, 97),
            DenseInterval::new(51, 81),
            DenseInterval::new(29, 87),
            DenseInterval::new(61, 92),
            DenseInterval::new(12, 70),
            DenseInterval::new(26, 54),
            DenseInterval::new(12, 52),
            DenseInterval::new(52, 54),
            DenseInterval::new(12, 72),
            DenseInterval::new(12, 58),
            DenseInterval::new(44, 88),
        ]);
        src.sort();
        let s = DenseStabby::new(100, &src);

        for qi in [
            DenseInterval::new(39, 84),
            DenseInterval::new(13, 44),
            DenseInterval::new(2, 57),
            DenseInterval::new(50, 75),
            DenseInterval::new(2, 11),
        ] {
            let mut expected: Vec<DenseInterval> = Vec::new();
            for ivl in src.iter() {
                if (ivl.first <= qi.last) && (ivl.last >= qi.first) {
                    expected.push(*ivl);
                }
            }
            println!("[{}, {}] -> {}", qi.first, qi.last, expected.len());
            assert_eq!(s.stab_interval(&qi), expected);
        }
    }

    #[test]
    fn test_stabby_4a() {
        let src: Vec<DenseInterval> = Vec::from([
            DenseInterval::new(0, 10),
            DenseInterval::new(0, 58),
            DenseInterval::new(2, 10),
            DenseInterval::new(2, 44),
            DenseInterval::new(4, 10),
            DenseInterval::new(4, 58),
            DenseInterval::new(6, 10),
            DenseInterval::new(6, 58),
            DenseInterval::new(8, 10),
            DenseInterval::new(8, 36),
            DenseInterval::new(12, 14),
            DenseInterval::new(16, 18),
            DenseInterval::new(20, 22),
            DenseInterval::new(24, 26),
            DenseInterval::new(28, 30),
            DenseInterval::new(32, 34),
            DenseInterval::new(32, 36),
            DenseInterval::new(38, 40),
            DenseInterval::new(42, 44),
            DenseInterval::new(46, 50),
            DenseInterval::new(48, 50),
            DenseInterval::new(52, 54),
            DenseInterval::new(56, 58),
            DenseInterval::new(60, 64),
            DenseInterval::new(60, 134),
            DenseInterval::new(60, 136),
            DenseInterval::new(60, 142),
            DenseInterval::new(62, 64),
            DenseInterval::new(62, 132),
            DenseInterval::new(66, 68),
            DenseInterval::new(70, 72),
            DenseInterval::new(74, 76),
            DenseInterval::new(78, 80),
            DenseInterval::new(82, 84),
            DenseInterval::new(86, 88),
            DenseInterval::new(90, 92),
            DenseInterval::new(94, 96),
            DenseInterval::new(98, 100),
            DenseInterval::new(102, 104),
            DenseInterval::new(106, 108),
            DenseInterval::new(110, 112),
            DenseInterval::new(114, 116),
            DenseInterval::new(118, 120),
            DenseInterval::new(122, 124),
            DenseInterval::new(126, 128),
            DenseInterval::new(130, 132),
            DenseInterval::new(130, 134),
            DenseInterval::new(130, 136),
            DenseInterval::new(138, 140),
            DenseInterval::new(144, 146),
            DenseInterval::new(144, 258),
            DenseInterval::new(144, 260),
            DenseInterval::new(148, 150),
            DenseInterval::new(152, 154),
            DenseInterval::new(156, 158),
            DenseInterval::new(160, 162),
            DenseInterval::new(164, 166),
            DenseInterval::new(168, 170),
            DenseInterval::new(172, 174),
            DenseInterval::new(176, 178),
            DenseInterval::new(180, 182),
            DenseInterval::new(184, 186),
            DenseInterval::new(188, 190),
            DenseInterval::new(192, 194),
            DenseInterval::new(196, 198),
            DenseInterval::new(200, 202),
            DenseInterval::new(204, 206),
            DenseInterval::new(208, 210),
            DenseInterval::new(212, 214),
            DenseInterval::new(216, 218),
            DenseInterval::new(220, 222),
            DenseInterval::new(224, 226),
            DenseInterval::new(228, 230),
            DenseInterval::new(232, 234),
            DenseInterval::new(236, 238),
            DenseInterval::new(240, 244),
            DenseInterval::new(242, 244),
            DenseInterval::new(242, 256),
            DenseInterval::new(246, 248),
            DenseInterval::new(250, 252),
            DenseInterval::new(254, 256),
            DenseInterval::new(254, 258),
        ]);

        let s = DenseStabby::new(261, &src);

        for q in 0..259 {
            let mut expected: Vec<DenseInterval> = Vec::new();
            for ivl in src.iter() {
                if ivl.first <= q && q <= ivl.last {
                    expected.push(*ivl);
                }
            }
            assert_eq!(s.stab(q), expected);
        }
    }

    #[test]
    fn test_stabby_5() {
        let src: Vec<DenseInterval> = Vec::from([
            DenseInterval::new(2, 4),
            DenseInterval::new(6, 8),
            DenseInterval::new(10, 12),
            DenseInterval::new(14, 16),
            DenseInterval::new(18, 20),
            DenseInterval::new(22, 24),
        ]);

        let s = DenseStabby::new(25, &src);

        assert_eq!(
            s.stab_interval(&DenseInterval {
                first: 13,
                last: 17
            }),
            vec![DenseInterval::new(14, 16)]
        );
    }
}
