extern crate clap;

use std::collections::HashMap;

use clap::{App, Arg, ArgMatches, SubCommand};

use crate::pwmlib::pwm::read_pwm_as_tab;

mod pwmlib;

fn main() {

    // Argument parser
    let matches: ArgMatches = App::new("pwm_phylo")
        .version("0.0.1")
        .author("Denis Puthier <puthier@gmail.com>")
        .about("Scan MAF files with PWM matrices.")
        .subcommand(SubCommand::with_name("score")
            .arg(Arg::with_name("converted-maf")
                .short("i")
                .long("converted-maf")
                .takes_value(true)
                .help("The MAF file converted using 'convert' command."))
            .arg(Arg::with_name("input-pwm")
                .short("p")
                .long("pwm")
                .takes_value(true)
                .help("The position weight matrix (not a count matrix !).").required(true))
            .arg(Arg::with_name("pwm-name")
                .short("n")
                .long("pwm-name")
                .takes_value(true)
                .help("The name of the pwm matrix.").required(true))
            .arg(Arg::with_name("p-value-threshold")
                .short("t")
                .long("p-value-threshold")
                .takes_value(true)
                .help("Motifs below with p-values threshold won't be printed.").required(true)))
        .subcommand(SubCommand::with_name("convert")
            .arg(Arg::with_name("input-maf")
                .short("i")
                .long("input-maf")
                .takes_value(true)
                .help("The input MAF file"))
            .arg(Arg::with_name("species-list")
                .short("l")
                .long("species-list")
                .takes_value(true)
                .help("A comma separated list of species.").required(true)))
        .get_matches();



    // You can handle information about subcommands by requesting their matches by name
    // (as below), requesting just the name used, or both at the same time
    if let Some(matches) = matches.subcommand_matches("score") {


        // The pwm as a hash -> (nucleotide, position) -> value

        let mut _pwm_to_score: HashMap<(char, i32), f64> = HashMap::new();

        // pwm size
        let _pwm_size: i32;


        // Path to the pwm
        let input_pwm_path;
        input_pwm_path = matches.value_of("input-pwm").unwrap().to_string();
        //input_pwm_path = matches.value_of("test").value_of("input-pwm").unwrap().to_string();
        let input_pwm_path_str = input_pwm_path.as_str();

        // P-value threshold
        let _score_treshold;
        _score_treshold = matches.value_of("p-value-threshold").unwrap().parse::<f64>().unwrap();


        // Name of the pwm
        let _pwm_name: String;
        _pwm_name = matches.value_of("pwm-name").unwrap().to_string();
        let _pwm_name_str = _pwm_name.as_str();

        // Read the PWM file
        // Return pwm_size
        // and fill _pwm_to_score

        _pwm_size = read_pwm_as_tab(input_pwm_path_str, &mut _pwm_to_score);


        pwmlib::maf::process_converted_maf(&_pwm_to_score,
                                         _pwm_size,
                                         _pwm_name_str,
                                         _score_treshold);



    } else if let Some(matches) = matches.subcommand_matches("convert") {


        // Get the list of species
        let species_list_str;
        species_list_str = matches.value_of("species-list").unwrap();
        let _species_list_vec: Vec<String> = species_list_str.split(",").map(|s| s.to_string()).collect();


        // Read the maf file.
        let _maf_file = matches.value_of("input-maf").unwrap_or("");

        if _maf_file.len() > 0 {
            pwmlib::maf::convert_maf(_maf_file,
                                        _species_list_vec.clone());
        } else {
            pwmlib::utils::message("MAF file could not be found.")
        }

    }


    }



