use std::collections::HashMap;

extern crate noodles;
extern crate stabby;

fn get_gene_name(rec: &noodles::gtf::Record) -> String {
    for e in rec.attributes().iter() {
        if e.key() == "gene_name" {
            return e.value().to_string();
        }
    }
    return "*".to_string();
}

fn main() -> std::io::Result<()> {
    let args = Vec::from_iter(std::env::args());
    if args.len() != 2 {
        println!("to run this example, get yourself a Gencode file. For example you could download one such as");
        println!("https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/gencode.v44.basic.annotation.gtf.gz");
        println!("then invoke the code in the following manner:");
        println!("    cargo run --example gencode -- gencode.v44.basic.annotation.gtf.gz");
    }
    let src = &args[1];
    let mut gtf = std::fs::File::open(src)
        .map(flate2::read::MultiGzDecoder::new)
        .map(std::io::BufReader::new)
        .map(noodles::gtf::Reader::new)?;

    let mut features = HashMap::new();
    let mut index = Vec::new();
    let mut l = 0;
    for res in gtf.records() {
        let rec = res?;
        let chrom = rec.reference_sequence_name();
        if chrom != "chr1" {
            continue;
        }
        if rec.end().get() > l {
            l = rec.end().get();
        }
        let ivl = stabby::Interval::new(rec.start().get() as u64, rec.end().get() as u64);
        let kind = rec.ty().to_string();
        let nm = get_gene_name(&rec);
        let feat_vec = features.entry(ivl).or_insert(Vec::new());
        feat_vec.push((nm, kind));
        index.push(ivl);
    }
    index.sort();
    index.dedup();
    let s = stabby::Stabby::new(&index);

    let mut depth = Vec::new();
    for q in 1..=l {
        let v = s.stab(q as u64);
        if depth.len() <= v.len() {
            depth.resize(v.len() + 1, 0);
        }
        depth[v.len()] += 1;
    }
    println!("{:?}", depth);
    Ok(())
}
