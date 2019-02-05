// compact representation of a kmer string
use bit_vec::BitVec;
use std::fmt;
use std::fmt::Debug;
use std::ops::Index;

#[derive(Clone)]
pub struct Kmer {
    contents: BitVec,
}

fn base2bibit(c : char) -> (bool, bool) {
    match c {
        'A' => (false, false),
        'T' => (false, true),
        'C' => (true, false),
        'G' => (true, true),
        _ => panic!("character not valid genome: {}", c),
    }
}

fn bibit2base(tp : (bool, bool)) -> char {
    match tp {
        (false, false) => 'A',
        (false, true) => 'T',
        (true, false) => 'C',
        (true, true) => 'G',
    }
}

impl Kmer {
    pub fn new() -> Kmer {
        Kmer {
            contents: BitVec::with_capacity(80 * 2),
        }
    }

    pub fn from_string(str: &str) -> Kmer {
        let mut v = BitVec::with_capacity(85 * 2);
        for c in str.chars() {
            let ele = base2bibit(c);
            v.push(ele.0);
            v.push(ele.1);
        }
        Kmer { contents: v }
    }

    pub fn push(&mut self, c: char) {
        let (l, r) = base2bibit(c);
        self.contents.push(l);
        self.contents.push(r);
    }

    pub fn rep(&self) -> String {
        let mut s = String::new();
        let n = self.contents.len() / 2;
        for i in 0..n {
            let c = bibit2base((self.contents[i * 2], self.contents[i * 2 + 1]));
            s.push(c);
        }
        s
    }

    pub fn is_edge(&self, rhs: &Kmer) -> bool {
        let tests = vec![5, 10, 15, 30, 40, 50, 60, 55];
        for t in tests {
            if self[t] != rhs[t - 1] {
                return false;
            }
        }
        return true;
    }
}

impl Index<usize> for Kmer {
    type Output = char;

    fn index(&self, i: usize) -> &char {
        match (self.contents[i * 2], self.contents[i * 2 + 1]) {
            (false, false) => &'A',
            (false, true) => &'T',
            (true, false) => &'C',
            (true, true) => &'G',
        }
    }
}
