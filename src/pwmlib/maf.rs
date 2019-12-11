extern crate regex;

use std::collections::HashMap;
use std::fs::File;
use std::io::BufRead;
use std::io::BufReader;
use std::io;
use regex::Regex;

pub fn convert_maf(input_maf: &str,
                   _species_list_vec: Vec<String>) {

    // declare several hashes
    let mut _species_to_seq: HashMap<String, String> = HashMap::new();
    let mut _species_to_chr: HashMap<String, String> = HashMap::new();
    let mut _species_to_strand: HashMap<String, String> = HashMap::new();
    let mut _species_to_start: HashMap<String, String> = HashMap::new();

    //The number of maf records processed
    let _nb_record = &mut 0;

    // This regex will be used to find
    // the sequences and extract, later,
    // all associated elements.
    let _re: Regex = Regex::new(r"s\s(\w+)\.(\w+)\s+(\d+)\s+(\d+)\s+([\+\-])\s+(\d+)\s+([AGCTagctnN\-]+)").unwrap();

    // Open the file
    self::super::utils::message(&format!("Found{} ", input_maf));
    let _file_handler = File::open(input_maf).unwrap();
    let _buffer_reader = BufReader::new(&_file_handler);

    // Loop through the _line.
    // Extract the required elements
    // (chrom, pos, species...).

    println!("##{}", _species_list_vec.join(","));

    for (_num, _result) in _buffer_reader.lines().enumerate() {
        let _line = _result.unwrap();

        if _line.starts_with("s ") {
            //'^s...' is a sequence.
            let caps = _re.captures(&_line).unwrap();
            let _species = caps.get(1).unwrap().as_str();
            let _chr = caps.get(2).unwrap().as_str();
            let _start = caps.get(3).unwrap().as_str();
            let _size = caps.get(4).unwrap().as_str();
            let _strand = caps.get(5).unwrap().as_str();
            let _src_size = caps.get(6).unwrap().as_str();
            let _seq = caps.get(7).unwrap().as_str();
            _species_to_seq.insert(_species.to_owned(), _seq.to_owned());
            _species_to_chr.insert(_species.to_owned(), _chr.to_owned());
            _species_to_strand.insert(_species.to_owned(), _strand.to_owned());
            _species_to_start.insert(_species.to_owned(), _start.to_owned());
        } else if _line.starts_with("a score=") {
            //'a score=' is the beginning of a record

            let mut to_process = 1;
            for s in _species_list_vec.clone() {
                if !_species_to_seq.contains_key(&s) {
                    to_process = 0;
                    break;
                }
            }

            if *_nb_record > 0 && to_process == 1 {

                for s in _species_list_vec.clone() {
                    print!("{}|", s);
                }
                for s in _species_list_vec.clone() {
                    print!("{}|", _species_to_chr[&s]);
                }
                for s in _species_list_vec.clone() {
                    print!("{}|", _species_to_start[&s]);
                }
                for s in _species_list_vec.clone() {
                    print!("{}|", _species_to_strand[&s]);
                }
                for s in _species_list_vec.clone() {
                    print!("{}|", _species_to_seq[&s]);
                }

                print!("\n");
            }
            _species_to_seq.clear();
            _species_to_chr.clear();
            _species_to_strand.clear();
            _species_to_start.clear();
            *_nb_record += 1;
        }
    }
}


//-----------------------------------------------------------------
//   Process the converted maf file
//-----------------------------------------------------------------

pub fn process_converted_maf(_pwm_to_score: &HashMap<(char, i32), f64>,
                               _pwm_size: i32,
                               _pwm_name: &str,
                               _score_treshold: f64) {

    // The number of motif/position tested
    let mut _counter:  i32 =  0;

    let stdin = io::stdin();
    let handle = stdin.lock();

    // Loop through the _line.
    // Extract the required elements
    // (chrom, pos, species...).


    let mut _species_list_vec: Vec<String> = std::vec::Vec::new();
    let mut _nb_species: i32 = 0;

    for (_num, _result) in handle.lines().enumerate() {
        let _line = _result.unwrap();



        // Retrieve the species list
        if _line.starts_with("##") {
            let _line = _line.replace("##", "");
            let _line = _line.replace("\n", "");
            let tmp:Vec<&str> = _line.split(",").collect();
            for  i in tmp{
                _species_list_vec.push(i.to_string());
            }
            _nb_species = _species_list_vec.len() as i32;

        }else{
            // Split the line.
            let _line_split: Vec<&str> = _line.split("|").collect();

            // Retrieve the sequence list
            let seq_idx = _nb_species * 4;
            let _seq_list = &_line_split[seq_idx as usize..];
            let seq_len = _seq_list[0].len() as i32;

            // Retrieve the chrom list
            let chr_idx_s = _nb_species as usize * 1;
            let chr_idx_e = chr_idx_s + _nb_species as usize ;
            let chr_list = &_line_split[chr_idx_s..chr_idx_e];

            // Retrieve the start list
            let start_idx_s = _nb_species as usize * 2;
            let start_idx_e = start_idx_s + _nb_species as usize ;
            let start_list = &_line_split[start_idx_s..start_idx_e];

            // Retrieve the strand list
            let strand_idx_s = _nb_species as usize * 3;
            let strand_idx_e = strand_idx_s + _nb_species as usize ;
            let strand_list = &_line_split[strand_idx_s..strand_idx_e];


            // The number of motif position
            let _nb_pwm_pos = seq_len - _pwm_size + 1;

            // Computing the motif score at position 'seq_pos'

            for _seq_start in 0.._nb_pwm_pos {

                let _seq_end = _seq_start + _pwm_size;

                let mut _motif_list:Vec<&str> = std::vec::Vec::new();
                let mut _motif_score:Vec<f64> = std::vec::Vec::new();
                let mut _motif_has_pos:Vec<i32> = std::vec::Vec::new();

                for _spe_pos in 0.._nb_species {

                    // Extract the motif at the given
                    // position
                    let _seq_start = _seq_start as usize;
                    let _seq_end = _seq_end as usize;
                    let _motif = &_seq_list[_spe_pos as usize][_seq_start.._seq_end];

                    // Score is set to 0
                    // at the moment
                    let mut score:f64 =0.0;

                    // Now iterate over motif chars
                    // to compute the associated score.
                    for (motif_pos, c) in _motif.chars().enumerate(){

                        // Switch char to lowercase
                        let c = c.to_lowercase();
                        let c = c.to_string().chars().next().unwrap();

                        // If there is an insertion/deletion
                        // score is set to -10000
                        if c == '-' || c == 'n' {
                            score = -100000.0;
                            break;
                        } else {

                            // Add the nucleotide score to the motif score
                            let tmp_pos = motif_pos as i32 ;
                            score += *_pwm_to_score.get(&(c, tmp_pos)).unwrap() ;
                        }


                    }

                    // Add the motif to the list species motifs
                    _motif_list.push(_motif);
                    _motif_score.push(score);
                    _counter += 1;

                    let my_string = format!("--| #Computed motifs: {}", _counter);
                    let nb_repeat = my_string.len();
                    eprint!("--| #Computed motifs: {}", _counter);
                    eprint!("{}", "\x08".to_string().repeat(nb_repeat));


                }
                // Write the result for the
                // corresponding motif.

                let mut output_line: Vec<String> = std::vec::Vec::new();

                // check whether score is
                // greater than threshold
                let mut mean_score :f64 = 0.0;
                for sc in _motif_score.clone() {
                    mean_score += sc;
                }
                mean_score = mean_score / _nb_species as f64;

                if mean_score < _score_treshold {
                    continue;
                }


                for _spe_pos in 0.._nb_species {
                    // Chromosome
                    output_line.push(chr_list[_spe_pos as usize].to_string());

                    // Start
                    let mut start_new = start_list[_spe_pos as usize].to_string().parse::<i32>().unwrap();
                    start_new += _seq_start;
                    let end_new = start_new + _pwm_size;
                    let start_new = format!("{}", start_new).to_string();
                    output_line.push(start_new);

                    // End
                    let end_new = format!("{}", end_new).to_string();
                    output_line.push(end_new);

                    // Name
                    let name_motif = format!("{}|{}",
                                             _species_list_vec[_spe_pos as usize].clone(),
                                             _motif_list[_spe_pos as usize]);
                    output_line.push(name_motif);

                    // Score
                    output_line.push(_motif_score[_spe_pos as usize].to_string());

                    // Strand
                    output_line.push(strand_list[_spe_pos as usize].to_string());


                }

                println!("{}", output_line.join("\t"));

            }


        }
    }
}

