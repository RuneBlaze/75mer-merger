// use bincode::{deserialize_from, serialize_into};
use itertools::join;
use serde::{Deserialize, Serialize};
use std::cmp::Ordering;
use std::collections::hash_map::DefaultHasher;
use hashbrown::HashMap;
use std::collections::VecDeque;
use std::error::Error;
use std::fs;
use std::fs::File;
use std::hash::{Hash, Hasher};
use std::io::BufRead;
use std::io::BufReader;
use std::io::Write;
mod kmer;
use std::env;
use self::kmer::Kmer;

// const DATASET: &str = "CC010_F4396";

type Path = Vec<u32>;
type Graph = HashMap<u32, Vec<u32>>;

fn calculate_hash<T: Hash>(t: &T) -> u64 {
    let mut s = DefaultHasher::new();
    t.hash(&mut s);
    s.finish()
}

fn dataset_name() -> String {
    return env::args().nth(1).unwrap();
}

#[derive(Debug, PartialEq, Eq, Clone, Serialize, Deserialize)]
struct Entry {
    id: u32,
    hash: u64,
    is_prefix: bool,
}

impl Ord for Entry {
    fn cmp(&self, other: &Entry) -> Ordering {
        (self.hash, self.is_prefix).cmp(&(other.hash, other.is_prefix))
    }
}

impl PartialOrd for Entry {
    fn partial_cmp(&self, other: &Entry) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

#[derive(Clone)]
struct Row {
    kmer: Kmer,
    pre_fcnt: i64,
    pre_rcnt: i64,
    fcnt: i64,
    rcnt: i64,
}

#[derive(Clone, Debug)]
struct Bucket {
    left: Vec<Entry>,
    right: Vec<Entry>,
}


impl Bucket {
    fn new() -> Bucket {
        Bucket {
            left: vec![],
            right: vec![],
        }
    }

    fn push(&mut self, entry: Entry) {
        if entry.is_prefix {
            self.right.push(entry);
        } else {
            self.left.push(entry);
        }
    }

    fn is_fit(&self) -> bool {
        !self.left.is_empty() && !self.right.is_empty()
    }

    fn is_empty(&self) -> bool {
        self.left.is_empty() && self.right.is_empty()
    }

    fn can_push(&self, entry: &Entry) -> bool {
        if self.is_empty() {
            return true;
        } else {
            if !self.left.is_empty() {
                return self.left[0].hash == entry.hash;
            }

            if !self.right.is_empty() {
                return self.right[0].hash == entry.hash;
            }

            return false;
        }
    }
}

fn mk_entries() -> Vec<Entry> {
    let mut res = vec![];
    let f = File::open(format!("{}.csv", dataset_name())).unwrap();
    let file = BufReader::new(&f);
    for (num, line) in file.lines().enumerate() {
        if num == 0 {
            continue;
        }
        let l = line.unwrap();
        let x = l.split(",").next().unwrap();
        let pre = &x[0..74];
        let suf = &x[1..75];
        let pre_h = calculate_hash(&pre);
        let suf_h = calculate_hash(&suf);
        res.push(Entry {
            id: (num - 1) as u32,
            hash: pre_h,
            is_prefix: true,
        });
        res.push(Entry {
            id: (num - 1) as u32,
            hash: suf_h,
            is_prefix: false,
        });
    }
    return res;
}

// TODO: this is an O(n lg n) algorithm.
// the thing is that, (unless we use a different hashmap)
// it seems that this O(n lg n) algorithm outperforms using hashmap for matching
fn equivalence_classes(coll: &mut Vec<Entry>) -> Vec<Bucket> {
    coll.sort_unstable();
    let mut buckets: Vec<Bucket> = vec![];
    let mut buf = Bucket::new();
    for ele in coll {
        if buf.can_push(ele) {
            buf.push(ele.clone());
        } else {
            if buf.is_fit() {
                buckets.push(buf.clone());
            }
            buf = Bucket::new();
        }
    }
    return buckets;
}

fn mk_graph(ecs: &Vec<Bucket>, rows: &Vec<Row>) -> HashMap<u32, Vec<u32>> {
    let mut res: HashMap<u32, Vec<u32>> = HashMap::new();
    for i in 0..(rows.len() as u32) {
        res.insert(i, vec![]);
    }
    for b in ecs {
        let lhs = &b.left;
        let rhs = &b.right;
        for l in lhs {
            for r in rhs {
                if is_true_edge(rows, l.id, r.id) {
                    res.get_mut(&l.id).unwrap().push(r.id);
                }
            }
        }
    }
    return res;
}

fn produce_cached_result() -> (Graph, Vec<Row>) {
    let mut entries = mk_entries();
    let ecs = equivalence_classes(&mut entries);
    let rows = load_rows();
    let graph = mk_graph(&ecs, &rows);
    // let serialized = serde_json::to_string(&graph).unwrap();
    // let serialized = serialize(&graph).unwrap();
    // let mut out = File::create(format!("graph_{}.bc", DATASET)).unwrap();
    // serialize_into(out, &graph).unwrap();
    return (graph, rows);
    // fs::write(format!("graph_{}.json", DATASET), serialized).expect("should output graph");
}

// fn load_graph() -> Result<Graph, Box<Error>> {
//     let file = File::open(format!("graph_{}.bc", DATASET))?;
//     let g: Graph = deserialize_from(file)?;
//     Ok(g)
// }

fn load_rows() -> Vec<Row> {
    let mut res = vec![];
    let f = File::open(format!("{}.csv", dataset_name())).unwrap();
    let file = BufReader::new(&f);
    for (num, line) in file.lines().enumerate() {
        if num == 0 {
            continue;
        }
        let l = line.unwrap();
        let mut iter = l.split(",");
        let x = Kmer::from_string(iter.next().unwrap());
        let n1 = iter.next().unwrap().parse::<i64>().unwrap();
        let n2 = iter.next().unwrap().parse::<i64>().unwrap();
        let n3 = iter.next().unwrap().parse::<i64>().unwrap();
        let n4 = iter.next().unwrap().parse::<i64>().unwrap();
        res.push(Row {
            kmer: x,
            pre_fcnt: n1,
            pre_rcnt: n2,
            fcnt: n3,
            rcnt: n4,
        });
    }
    return res;
}

fn is_true_edge(rows: &Vec<Row>, l: u32, r: u32) -> bool {
    rows[(l) as usize].kmer.is_edge(&rows[(r) as usize].kmer)
}

fn dfs(rows: &Vec<Row>, g: &Graph) -> Vec<Path> {
    let mut res: Vec<Path> = vec![];
    let mut stack: Vec<u32> = vec![];
    let mut buf: Vec<u32> = vec![];
    let mut parent: HashMap<u32, u32> = HashMap::new();
    let mut indegree: HashMap<u32, i32> = HashMap::new();
    for i in 0..(rows.len() as u32) {
        indegree.insert(i, 0);
    }
    for (_k, v) in g.iter() {
        // indegree.entry(*k).or_insert(0);
        for target in v {
            let e = indegree.entry(*target).or_insert(0);
            *e += 1;
        }
    }

    for (k, v) in indegree.iter() {
        if *v == 0 {
            stack.push(*k);
        }
    }

    // println!("{:?}", indegree);

    let def = vec![];
    while let Some(cur) = stack.pop() {
        let next = g.get(&cur).unwrap_or(&def);
        if next.is_empty() {
            // the path ends here
            let mut pt = cur;
            let mut stk = vec![];
            while let Some(prev) = parent.get(&pt) {
                stk.push(pt);
                pt = *prev;
            }
            stk.push(pt);
            while let Some(e) = stk.pop() {
                buf.push(e);
            }
            res.push(buf.clone());
            buf.clear();
        } else {
            for n in next.iter() {
                stack.push(*n);
                parent.insert(*n, cur);
            }
        }
    }

    return res;
}

impl Row {
    fn rep(&self) -> String {
        format!(
            "{},{},{},{},{}",
            self.kmer.rep(),
            self.pre_fcnt,
            self.pre_rcnt,
            self.fcnt,
            self.rcnt
        )
    }
}

fn path_to_string(rows: &Vec<Row>, path: &Path) -> String {
    let mut it = path.iter();
    let mut entry = rows[(*it.next().unwrap()) as usize].clone();
    while let Some(n) = it.next() {
        let r = &rows[(*n) as usize];
        entry.kmer.push(r.kmer[74]);
    }
    let mut r = entry.rep();
    r.push_str(&format!(",{}", join(path, "-")));
    r
}

fn process_output(rows: &Vec<Row>, paths: &Vec<Path>) -> () {
    let mut output = File::create(format!("{}.merged.csv", dataset_name())).unwrap();
    writeln!(output, "sequence,lo,hi,pre45cnt,count,path").unwrap();
    for p in paths {
        let s = path_to_string(rows, p);
        writeln!(output, "{}", s).unwrap();
    }
}

fn main() {
    // produce_cached_result();
    // let mut entries = mk_entries();
    // let ecs = equivalence_classes(&mut entries);
    // let graph = mk_graph(&ecs);
    // for (k, v) in graph.iter() {
    //     println!("{};{}",k,join(v, ";"));
    // }
    // let serialized = serialize(&graph).unwrap();
    // let mut out = File::create("graph_WEB_FNNN.bc").unwrap();
    // serialize_into(out, &graph).unwrap();
    // fs::write("graph_WEB_FNNN.json", serialized).expect("should output graph");
    if env::args().nth(1).is_none() {
        panic!("please pass the dataset name as the first argument, e.g. <dataset_name>.csv must exist.");
    }

    println!("process started");
    let (graph, rows) = produce_cached_result();
    println!("producing cached result");
    // let graph = load_graph().unwrap();
    // let rows = load_rows();
    let mut paths = dfs(&rows, &graph);
    println!("sorting");
    paths.sort_unstable();
    // println!("# of paths {}", paths.len());
    println!("sorting finished");

    let mut cnt: i64 = 0;
    for p in paths.iter() {
        cnt += p.len() as i64;
    }
    println!("{}", cnt);
    process_output(&rows, &paths);
    println!("# of rows, previously: {}", rows.len());
    println!("# of rows, after: {}", paths.len());
    println!("after / prev: {}", paths.len() as f64 / rows.len() as f64);
    let maxlen = paths.iter().max_by(|x, y| x.len().cmp(&y.len())).unwrap().len();
    println!("max length of paths: {}", maxlen);
    // for i in 0..10 {
    //     for p in &paths[i] {
    //         println!("{:?}", rows[*p as usize].kmer.rep());
    //     }
    //     println!("========");
    // }

    // let mut cnt = 0;
    // for p in paths.iter() {
    //     if p.len() == 1 {
    //         cnt += 1;
    //         // println!("{:?}", p);
    //     }
    // }

    // println!("{}", cnt);

    // let mut m = Kmer::from_string("ACTGTCG");
    // println!("{}", m.rep());
    // m.push('C');
    // println!("{}", m.rep());
    // let mut eg_graph : HashMap<u32, Vec<u32>> = HashMap::new();
    // eg_graph.insert(0, vec![1]);
    // eg_graph.insert(1, vec![3]);
    // eg_graph.insert(3, vec![2,7]);
    // eg_graph.insert(5, vec![4]);
    // eg_graph.insert(8, vec![4,6]);
    // eg_graph.insert(6, vec![7]);

    // let r = dfs(&eg_graph);
    // println!("Graph Reconstructed");
    // let rows = load_rows();
    // println!("Rows Loaded");
    // let mut cnt = 0;
    // for (k, v) in g.iter() {
    //     for target in v {
    //         let r = rows[(*k) as usize].kmer.is_edge(&rows[(*target) as usize].kmer);
    //         if !r {
    //             println!("{} -> {}", k, target);
    //             cnt+=1;
    //         }
    //         println!("{}", rows[(*k) as usize].kmer.rep());
    //         println!("{}", &rows[(*target) as usize].kmer.rep());
    //         break;
    //     }
    //     break;
    // }
    // println!("total mishashed: {}", cnt);
    // println!("{:?}", r);
}
