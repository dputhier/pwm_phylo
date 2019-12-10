extern crate regex;

use std::collections::HashMap;
use std::fs::File;
use std::io::BufRead;
use std::io::BufReader;

use regex::Regex;

//-----------------------------------------------------------------
//   Process the maf file
//-----------------------------------------------------------------

pub fn process_maf(input_maf: &str,
                   _pwm_to_score: HashMap<(char, i32), f64>,
                   _pwm_size: i32,
                   _pwm_name: &str,
                   _species_list_vec: Vec<String>,
                   _score_treshold: f64) {

    // declare several hashes
    let mut _species_to_seq: HashMap<String, String> = HashMap::new();
    let mut _species_to_chr: HashMap<String, String> = HashMap::new();
    let mut _species_to_strand: HashMap<String, String> = HashMap::new();
    let mut _species_to_start: HashMap<String, String> = HashMap::new();

    // The number of motif/position tested
    let _counter: &mut i32 = &mut 0;

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

    for (_num, _result) in _buffer_reader.lines().enumerate() {
        let _line = _result.unwrap();

        process_one_line(_line,
                         &_re,
                         &mut _species_to_seq,
                         &mut _species_to_chr,
                         &mut _species_to_strand,
                         &mut _species_to_start,
                         &_species_list_vec,
                         &_pwm_to_score,
                         _pwm_size,
                         _pwm_name,
                         _nb_record,
                         _counter,
                         _score_treshold);
    }
}


//-----------------------------------------------------------------
//   Extract the required fields / _line
//   Put them into corresponding hashes
//   Call the function that scan PWM
//-----------------------------------------------------------------

fn process_one_line(_line: String,
                    _re: &Regex,
                    _species_to_seq: &mut HashMap<String, String>,
                    _species_to_chr: &mut HashMap<String, String>,
                    _species_to_strand: &mut HashMap<String, String>,
                    _species_to_start: &mut HashMap<String, String>,
                    _species_list_vec: &Vec<String>,
                    _pwm_to_score: &HashMap<(char, i32), f64>,
                    _pwm_size: i32,
                    _pwm_name: &str,
                    _nb_record: &mut i32,
                    _counter: &mut i32,
                    _score_treshold: f64) {
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

        let mut to_process = 0;
        for s in _species_list_vec.clone() {
            if !_species_to_seq.contains_key(&s) {
                to_process = 1;
            }
        }
        if *_nb_record > 0 && to_process == 0 {
            scan_seqs(_species_to_seq,
                      _species_to_chr,
                      _species_to_strand,
                      _species_to_start,
                      _pwm_to_score,
                      _pwm_size,
                      _pwm_name,
                      _species_list_vec,
                      _score_treshold,
                      _counter);
        }
        _species_to_seq.clear();
        *_nb_record += 1;
    }
}


//-----------------------------------------------------------------
//   Scan the PWM over a record (a set of _lines 's=...')
//-----------------------------------------------------------------

fn scan_seqs(_species_to_seq: &mut HashMap<String, String>,
             _species_to_chr: &mut HashMap<String, String>,
             _species_to_strand: &mut HashMap<String, String>,
             _species_to_start: &mut HashMap<String, String>,
             __pwm_to_score: &HashMap<(char, i32), f64>,
             __pwm_size: i32,
             __pwm_name: &str,
             _species_list_vec: &Vec<String>,
             _score_treshold: f64,
             _counter: &mut i32) {

// species -> (chr, motif, pos, score)
    let mut species_to_motif: HashMap<String, Vec<(String, i32, String, f64, i8)>> = HashMap::new();

// the number of matrix position in the seq
    let mut nb_pwm_pos = 0;
//iterating over species
    for (_species, _seq) in _species_to_seq.iter() {
        let seq_len = _seq.len() as i32;
        nb_pwm_pos = seq_len - __pwm_size + 1;


// iterating over the sequence
        for seq_pos in 0..nb_pwm_pos {


            //to store the corresponding motif
            let mut _motif = "".to_string();
            //to store the corresponding motif score
            let mut _cur_species_pwm_score: f64 = 0.0;

            let mut has_indel: i8 = 0;
            eprint!("--| #Computed motifs: {}", *_counter);

            let my_string = format!("--| #Computed motifs: {}", *_counter);

            let nb_repeat = my_string.len();
            eprint!("{}", "\x08".to_string().repeat(nb_repeat));

            // Computing the motif score at position 'seq_pos'
            for pwm_pos in 0..__pwm_size - 1 {

                // the shift regarding seq_pos is 'pwm_pos'
                let nuc_pos = seq_pos as usize + pwm_pos as usize;

                // get the corresponding character
                let seq_char_raw = _seq.chars().nth(nuc_pos).unwrap();

                // add the char to the motif
                _motif.push_str(&seq_char_raw.to_string());

                // Convert char to lowercase
                let seq_char_to_low = seq_char_raw.to_lowercase();
                let seq_char = seq_char_to_low.to_string().chars().next().unwrap();


                //Check the char is not '-' or 'n' (N lowercased)
                if seq_char == '-' || seq_char == 'n' {
                    _cur_species_pwm_score = 0.0;
                    has_indel = 1;
                } else {
                    // Add the nucleotide score to the motif score
                    _cur_species_pwm_score += *__pwm_to_score.get(&(seq_char, pwm_pos + 1)).unwrap();
                }
            }

            *_counter += 1;

            species_to_motif.entry(_species.to_string()
            ).or_insert(Vec::new()
            ).push((
                _species_to_chr[&_species.to_owned()].clone(),
                _species_to_start[&_species.to_owned()].parse::<i32>().unwrap() + seq_pos,
                _motif.clone(),
                _cur_species_pwm_score,
                has_indel))
        }
    }

    let mut pos: i32 = 0;

    let mut nb_species = 0.0;

    for _ in _species_list_vec.clone() {
        nb_species += 1.0;
    }

    while pos < nb_pwm_pos {
        let mut has_indel: i8 = 0;
        let mut output_list: Vec<String> = Vec::new();
        let mut sum_pwm_scores: f64 = 0.0;
        for s in _species_list_vec.clone() {
            let v = &species_to_motif[&s][pos as usize];

            sum_pwm_scores += v.3;

            output_list.extend(vec![s,
                                    v.0.clone(),
                                    v.1.to_string(),
                                    v.2.clone(), format!("{:.4}", v.3.to_string())].iter().cloned());
            if v.4 == 1 {
                has_indel = 1;
            }
        }

        if sum_pwm_scores / nb_species > _score_treshold && has_indel == 0 {
            let output_str = output_list.join("\t");
            println!("{}", output_str);
        }
        pos += 1;
    }
}

