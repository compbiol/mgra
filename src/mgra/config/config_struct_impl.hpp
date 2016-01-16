#ifndef CONFIG_STRUCT_IMPL_HPP
#define CONFIG_STRUCT_IMPL_HPP

namespace config {

/*
 * Read configure file in JSON format and init all config
 */
template<class mcolor_t>
void load(mgra_config <mcolor_t> &cfg, std::istream &source) {
    INFO("Start load cfg file in MGRA format");

    Json::Value root;
    Json::Reader reader;
    bool isError = reader.parse(source, root);


    if (!isError) {
        std::string error_msg = "Failed to parse JSON config\n" + reader.getFormattedErrorMessages();
        std::cerr << error_msg << std::endl;
        exit(1);
    }

    cfg.load(root);

    INFO("Finish load cfg file in MGRA format")
}

/*
 * Convert string to color. If string is bad, terminate program.
 */
template<class mcolor_t>
mcolor_t mgra_config<mcolor_t>::name_to_mcolor(std::string const &temp) const {
    mcolor_t answer;
    if (temp[0] == '{' && temp[temp.length() - 1] == '}') {
        std::string current = "";
        for (size_t i = 1; i < temp.length(); ++i) {
            if (temp[i] == ',' || temp[i] == '}') {
                if (genome_number.find(current) != genome_number.end()) {
                    answer.insert(genome_number.find(current)->second);
                }
                current = "";
            } else if (std::isalpha(temp[i]) || std::isdigit(temp[i])) {
                current += temp[i];
            } else {
                ERROR("Bad format target " << temp);
                answer.clear();
                break;
            }
        }
    } else {
        ERROR("Bad format target " << temp);
    }
    return answer;
}

/*
 * Convert multicolor to string with genomes name.
 */
template<class mcolor_t>
std::string mgra_config<mcolor_t>::mcolor_to_name(mcolor_t const &color) const {
    /*if (multicolor_to_name.find(color) != multicolor_to_name.end()) {
        return multicolor_to_name.find(color)->second;
    } else {*/
    if (color.size() == 1) {
        return priority_name[color.cbegin()->first];
    }
    std::string answer = "{";
    if (!color.empty()) {
        std::string const &main_sym = priority_name[color.cbegin()->first];
        answer += main_sym;
        for (size_t i = 1; i < color.cbegin()->second; ++i) {
            answer += ("," + main_sym);
        }

        for (auto col = (++color.cbegin()); col != color.cend(); ++col) {
            std::string const &sym = priority_name[col->first];
            for (size_t i = 0; i < col->second; ++i) {
                answer += ("," + sym);
            }
        }
    }
    answer += '}';
    return answer;
    //}
}

/*
 * Different function, which load input.
 */
template<class mcolor_t>
void mgra_config<mcolor_t>::load(Json::Value const &root) {
    if (!root["genomes"].isNull()) {
        load_genomes(root["genomes"]);
    } else {
        std::cerr << "ERROR: genomes section is required" << std::endl;
        exit(1);
    }

    if (!root["format_of_blocks"].isNull()) {
        load_block_type(root["format_of_blocks"]);
    } else {
        std::cerr << "ERROR: format_of_blocks section is required" << std::endl;
        exit(1);
    }

    if (!root["files"].isNull()) {
        load_files(root["files"]);
    } else {
        exit(1);
    }

    if (!root["tree"].isNull()) {
        load_tree(root["tree"]);
    } else {
        std::cerr << "ERROR: tree section is required" << std::endl;
        exit(1);
    }

    if (!root["wgd"].isNull()) {
        load_wgd_events(root["wgd"]);
        structure::phyl_tree::refine_tree_by_wgds(phylotree, wgds_events);
    }

    if (!root["target"].isNull()) {
        load_target(root["target"]);
    }

    if (!root["Algorithm"].isNull()) {
        load_algorithm(root["Algorithm"]);
    } else if (how_build == default_algo) {
        default_algorithm();
    } else if (how_build == target_algo) {
        default_target_algorithm();
    } else if (how_build == wgd_algo) {
        default_wgd_algorithm();
    }

    if (!out_path_directory.empty()) {
        default_directory_organization();
    } else if (!root["output_directory"].isNull()) {
        load_output_directory(root["output_directory"]);
        default_directory_organization();
    } else {
        std::cerr << "ERROR: output directory is required. Please type -o <dirname> or in config file" << std::endl;
        exit(1);
    }

    if (!root["complections"].isNull()) {
        load_complections(root["target"]);
    }

    if (!root["saves_enable"].isNull()) {
        load_saves(root["saves_enable"]);
    }

    if (!root["debug_enable"].isNull()) {
        load_debug(root["debug_enable"]);
    }

    default_rgb_colors();
}

template<class mcolor_t>
Json::Value mgra_config<mcolor_t>::save() const {
    Json::Value root;

    root["genomes"] = save_genomes();

    root["format_of_blocks"] = save_block_type();

    root["files"] = save_files();

    root["tree"] = save_tree();

    if (!wgds_events.empty()) {
        root["wgd"] = save_wgd_events();
    }

    if (!target.empty()) {
        root["target"] = save_target();
    }

    root["Algorithm"] = save_algorithm();

    root["output_directory"] = save_output_directory();

    if (!completion.empty()) {
        root["complections"] = save_complections();
    }

    root["saves_enable"] = save_saves();

    root["debug_enable"] = save_debug();

    return root;
}

/**
 * Load functions
 */
template<class mcolor_t>
void mgra_config<mcolor_t>::load_algorithm(Json::Value const &algo) {
    if (!algo.isObject()) {
        std::cerr << "ERROR: algorithm must have object type" << std::endl;
    }

    // Parse rounds section
    if (!algo["rounds"].isNull()) {
        if (!algo["rounds"].isInt()) {
            std::cerr << "ERROR: rounds field in algorithm section must have int type" << std::endl;
            exit(1);
        }

        if (algo["rounds"].asInt() < 1 || algo["rounds"].asInt() > 3) {
            std::cerr << "ERROR: rounds field in algorithm section must have range between 1 and 3" << std::endl;
            exit(1);
        }

        rounds = (size_t) algo["rounds"].asInt();
    } else {
        std::cerr << "ERROR: algorithm section must have rounds field" << std::endl;
        exit(1);
    }

    // Parse pipeline section
    if (!algo["stages"].isNull()) {
        load_pipeline(algo["stages"]);
    } else {
        std::cerr << "ERROR: algorithm section must have stages field" << std::endl;
        exit(1);
    }
}

template<class mcolor_t>
void mgra_config<mcolor_t>::load_pipeline(Json::Value const &stages) {
    if (!stages.isArray()) {
        std::cerr << "ERROR: stages field in algorithm section must have array type" << std::endl;
        exit(1);
    }

    if (!stages.empty()) {
        std::cerr << "ERROR: stages field in algorithm section is array with at least one stage" << std::endl;
        exit(1);
    }

    for (auto const &stage : stages) {
        if (stage.asString() == "balance") {
            pipeline.push_back(kind_stage::balance_k);
        } else if (stage.asString() == "simple_path") {
            pipeline.push_back(kind_stage::simple_path_k);
        } else if (stage.asString() == "four_cycles") {
            pipeline.push_back(kind_stage::four_cycles_k);
        } else if (stage.asString() == "fair_edge") {
            pipeline.push_back(kind_stage::fair_edge_k);
        } else if (stage.asString() == "clone") {
            pipeline.push_back(kind_stage::clone_k);
        } else if (stage.asString() == "fair_clone_edge") {
            pipeline.push_back(kind_stage::fair_clone_edge_k);
        } else if (stage.asString() == "components") {
            pipeline.push_back(kind_stage::components_k);
        } else if (stage.asString() == "bruteforce") {
            if (how_build == target_algo) {
                std::cerr << "ERROR: Don't use bruteforce stage in target reconstruction" << std::endl;
                exit(1);
            }
            pipeline.push_back(kind_stage::bruteforce_k);
        } else if (stage.asString() == "blossomv") {
            if (how_build == target_algo) {
                std::cerr << "ERROR : Don't use blossom V stage in target reconstruction" << std::endl;
                exit(1);
            }
            pipeline.push_back(kind_stage::blossomv_k);
        } else {
            std::cerr << "Unknown option " << stage.asString() << std::endl;
            exit(1);
        }
    }

    size_component_in_bruteforce = 20;
}

template<class mcolor_t>
void mgra_config<mcolor_t>::load_wgd_events(Json::Value const &wgds) {
    if (!wgds.isArray()) {
        std::cerr << "ERROR: wgd section is array at least one wgd event" << std::endl;
        exit(1);
    }

    for (auto const &wgd : wgds) {
        load_wgd_event(wgd);
    }
}

template<class mcolor_t>
void mgra_config<mcolor_t>::load_wgd_event(Json::Value const &wgd) {
    if (!wgd.isObject()) {
        std::cerr << "ERROR: wgd must have object type" << std::endl;
        exit(1);
    }

    // Parse parent of wgd event
    if (wgd["parent"].isNull()) {
        std::cerr << "ERROR: wgd branch must have parent name" << std::endl;
        exit(1);
    }

    if (!wgd["parent"].isString()) {
        std::cerr << "ERROR: parent name in wgd branch must have string type." << std::endl;
        exit(1);
    }

    // Parse children of wgd event
    if (wgd["children"].isNull()) {
        std::cerr << "ERROR: wgd branch must have children name" << std::endl;
        exit(1);
    }

    if (!wgd["children"].isString()) {
        std::cerr << "ERROR: children name in wgd branch must have string type. " << std::endl;
        exit(1);
    }

    std::vector<std::string> names;
    if (!wgd["names"].isNull()) {
        if (!wgd["names"].isArray()) {
            std::cerr << "ERROR: wgs's names must have array type." << std::endl;
            exit(1);
        }

        if (wgd["names"].size() < 1) {
            std::cerr << "ERROR: wgd branch must have at least one name" << std::endl;
            exit(1);
        }

        for (auto name : wgd["names"]) {
            if (!name.isString()) {
                std::cerr << "ERROR: wgd's names must have string type" << std::endl;
                exit(1);
            }

            names.push_back(name.asString());
        }
    } else {
        std::cerr << "ERROR: wgd branch must have at least one name" << std::endl;
        exit(1);
    }

    wgds_events.insert(std::make_pair(std::make_pair(wgd["parent"].asString(), wgd["children"].asString()),
                       wgd_t(wgd["parent"].asString(), wgd["children"].asString(), names)));
}

template<class mcolor_t>
void mgra_config<mcolor_t>::load_target(Json::Value const &target_json) {
    if (!target_json.isString()) {
        std::cerr << "ERROR: target must have string type" << std::endl;
        exit(1);
    }

    target = target_json.asString();
    boost::trim(target);
}

template<class mcolor_t>
void mgra_config<mcolor_t>::load_complections(Json::Value const &twobreaks) {
    if (!twobreaks.isArray()) {
        std::cerr << "ERROR: complection section is array at least one two-break event" << std::endl;
        exit(1);
    }

    for (auto const &twobreak : twobreaks) {
        load_complection(twobreak);
    }
}

template<class mcolor_t>
void mgra_config<mcolor_t>::load_complection(Json::Value const &twobreak) {
    if (!twobreak.isString()) {
        std::cerr << "ERROR: twobreak must have string type" << std::endl;
        exit(1);
    }

    std::vector<std::string> mc(5);
    std::istringstream is(twobreak.asString());
    is >> mc[0] >> mc[1] >> mc[2] >> mc[3] >> mc[4];
    completion.push_back(twobreak_t(mc[0], mc[1], mc[2], mc[3], name_to_mcolor(mc[4])));
}

template<class mcolor_t>
void mgra_config<mcolor_t>::load_saves(Json::Value const &enable_saves) {
    if (!enable_saves.isBool()) {
        std::cerr << "ERROR: saves field must have bool type" << std::endl;
        exit(1);
    }

    is_saves = enable_saves.asBool();
}

template<class mcolor_t>
void mgra_config<mcolor_t>::load_debug(Json::Value const &enable_debug) {
    if (!enable_debug.isBool()) {
        std::cerr << "ERROR: debug field must have bool type" << std::endl;
        exit(1);
    }

    is_debug = enable_debug.asBool();
}

template<class mcolor_t>
Json::Value mgra_config<mcolor_t>::save_algorithm() const {
    Json::Value algo(Json::objectValue);
    algo["rounds"] = Json::Value((int) rounds);
    algo["stages"] = save_pipeline();
    return algo;
}

template<class mcolor_t>
Json::Value mgra_config<mcolor_t>::save_pipeline() const {
    Json::Value stages(Json::arrayValue);

    for (auto const &stage : pipeline) {
        if (stage == kind_stage::balance_k) {
            stages.append(Json::Value("balance"));
        } else if (stage == kind_stage::simple_path_k) {
            stages.append(Json::Value("simple_path"));
        } else if (stage == kind_stage::four_cycles_k) {
            stages.append(Json::Value("four_cycles"));
        } else if (stage == kind_stage::fair_edge_k) {
            stages.append(Json::Value("fair_edge"));
        } else if (stage == kind_stage::clone_k) {
            stages.append(Json::Value("clone"));
        } else if (stage == kind_stage::fair_clone_edge_k) {
            stages.append(Json::Value("fair_clone_edge"));
        } else if (stage == kind_stage::components_k) {
            stages.append(Json::Value("components"));
        } else if (stage == kind_stage::bruteforce_k) {
            stages.append(Json::Value("bruteforce"));
        } else if (stage == kind_stage::blossomv_k) {
            stages.append(Json::Value("blossomv"));
        }
    }

    return stages;
}

template<class mcolor_t>
Json::Value mgra_config<mcolor_t>::save_wgd_events() const {
    Json::Value events(Json::arrayValue);
    //TODO here
    //for (size_t ind = 0; ind < completion.size(); ++ind) {
    //    ;//events.append(save_wgd_event(ind));
    //}

    return events;
}

template<class mcolor_t>
Json::Value mgra_config<mcolor_t>::save_wgd_event(size_t ind) const {
    assert(ind != 0); //TODO here
    Json::Value wgd;
    /*wgd["parent"] = Json::Value(wgds_events[ind].get_parent());
    wgd["children"] = Json::Value(wgds_events[ind].get_children());
    wgd["names"] = Json::Value(Json::arrayValue);

    for (auto const &name : wgds_events[ind]) {
        wgd["names"].append(Json::Value(name));
    }*/

    return wgd;
}

template<class mcolor_t>
Json::Value mgra_config<mcolor_t>::save_target() const {
    return Json::Value(target);
}

template<class mcolor_t>
Json::Value mgra_config<mcolor_t>::save_complections() const {
    Json::Value complects(Json::arrayValue);

    for (auto it = completion.cbegin(); it != completion.cend(); ++it) {
        complects.append(save_complection(it));
    }

    return complects;
}

template<class mcolor_t>
Json::Value mgra_config<mcolor_t>::save_complection(
        typename std::list<twobreak_t>::const_iterator const &twobreak) const {
    std::ostringstream out;
    out << twobreak->get_vertex(0) << " " << twobreak->get_vertex(1) << " " << twobreak->get_vertex(2) << " " <<
    twobreak->get_vertex(3) << " " << mcolor_to_name(twobreak->get_multicolor());
    return Json::Value(out.str());
}

template<class mcolor_t>
Json::Value mgra_config<mcolor_t>::save_saves() const {
    return Json::Value(is_saves);
}

template<class mcolor_t>
Json::Value mgra_config<mcolor_t>::save_debug() const {
    return Json::Value(is_debug);
}

/**
 * Function, which initilization structure on default parameters
 */
template<class mcolor_t>
void mgra_config<mcolor_t>::default_algorithm() {
    size_component_in_bruteforce = 0;
    rounds = 3;
    pipeline.push_back(kind_stage::balance_k);
    pipeline.push_back(kind_stage::simple_path_k);
    pipeline.push_back(kind_stage::four_cycles_k);
    pipeline.push_back(kind_stage::fair_edge_k);
    pipeline.push_back(kind_stage::clone_k);
    pipeline.push_back(kind_stage::components_k);
    pipeline.push_back(kind_stage::change_canform_k);
    pipeline.push_back(kind_stage::blossomv_k);
}

template<class mcolor_t>
void mgra_config<mcolor_t>::default_wgd_algorithm() {
    size_component_in_bruteforce = 0;
    rounds = 1;
    pipeline.push_back(kind_stage::simple_path_k);
    pipeline.push_back(kind_stage::irregular_fair_edge_k);
}

template<class mcolor_t>
void mgra_config<mcolor_t>::default_target_algorithm() {
    size_component_in_bruteforce = 0;
    rounds = 3;
    pipeline.push_back(kind_stage::balance_k);
    pipeline.push_back(kind_stage::simple_path_k);
    pipeline.push_back(kind_stage::four_cycles_k);
    pipeline.push_back(kind_stage::fair_edge_k);
    pipeline.push_back(kind_stage::clone_k);
    pipeline.push_back(kind_stage::components_k);
    pipeline.push_back(kind_stage::change_canform_k);
}

template<class mcolor_t>
void mgra_config<mcolor_t>::default_directory_organization() {
    std::string const LOGGER_FILENAME = "mgra.log";
    std::string const INPUT_DIRNAME = "input";
    std::string const DEBUG_DIRNAME = "debug";
    std::string const SAVES_DIRNAME = "saves";
    std::string const GENOMES_DIRNAME = "genomes";
    std::string const TRANS_DIRNAME = "transformations";

    out_path_to_logger_file = path::append_path(out_path_directory, LOGGER_FILENAME);
    out_path_to_input_dir = path::append_path(out_path_directory, INPUT_DIRNAME);
    out_path_to_debug_dir = path::append_path(out_path_directory, DEBUG_DIRNAME);
    out_path_to_saves_dir = path::append_path(out_path_directory, SAVES_DIRNAME);
    out_path_to_genomes_dir = path::append_path(out_path_directory, GENOMES_DIRNAME);
    out_path_to_transfomations_dir = path::append_path(out_path_directory, TRANS_DIRNAME);
}

template<class mcolor_t>
void mgra_config<mcolor_t>::default_rgb_colors() {
    if (priority_name.size() < 10) {
        colorscheme = "set19";
        for (size_t i = 1; i < priority_name.size() + 1; ++i) {
            RGBcolors.push_back(std::to_string(i));
        }
        RGBcoeff = 1;
    } else if (priority_name.size() < 13) {
        colorscheme = "";
        size_t const number_colors = 13;
        std::string cols[number_colors] = {
                "red3", "green3", "blue", "purple", "black", "orange", "greenyellow", "pink",
                "cyan", "magenta", "yellow3", "grey70", "darkorange4"
        };

        for (size_t i = 0; i < number_colors; ++i) {
            RGBcolors.push_back(cols[i]);
        }

        RGBcoeff = 1;
    } else {
        colorscheme = "";
        size_t const number_colors = 136;
        std::string cols[number_colors] = {
                "#C91F16", "#CA2316", "#CC2A13", "#D03815", "#CC3615", "#D23D16", "#D34810", "#D44B0D",
                "#D9580E", "#D95B0F", "#DA5F0E", "#DB640D", "#DD680B", "#DC6E0D", "#E37509", "#E7860B",
                "#E68601", "#ED9E0A", "#F1A60A", "#F1AD06", "#EEAD00", "#F7BD00", "#F5B600", "#F9C60F",
                "#FCCF03", "#FBD50A", "#FCE503", "#F6E60A", "#EDE400", "#E6E209", "#DFDC09", "#D8DB09",
                "#CED500", "#C5D300", "#BED20F", "#B5CC0A", "#B0C50F", "#A5C50F", "#9FC40D", "#98C100",
                "#8FBB0C", "#92C00E", "#7FB513", "#80B812", "#77B312", "#71B30E", "#6FB20F", "#69AC12",
                "#69B011", "#55A41C", "#5AAA1D", "#4EA41D", "#47A41E", "#3A9B20", "#349B26", "#339C1E",
                "#2D9C1F", "#2B9B22", "#209426", "#1B9324", "#189425", "#039331", "#008C33", "#098D3A",
                "#088343", "#0C8C4B", "#00844A", "#048D5C", "#0E8C62", "#088D6C", "#008C7C", "#089394",
                "#05949D", "#0B8D94", "#0994A6", "#009DBE", "#0894B6", "#0D9BC4", "#009DCF", "#159BD4",
                "#1F9AD7", "#148DC6", "#1D93D1", "#147CB5", "#0D84C4", "#0C7BBB", "#0D73B3", "#106CAC",
                "#0363A3", "#135A9E", "#085293", "#0A5397", "#094A91", "#0D4284", "#06438A", "#0D3981",
                "#113279", "#123179", "#152C74", "#182A71", "#182369", "#1B1A5F", "#181C67", "#1D1A62",
                "#1D1259", "#1E1059", "#230E59", "#300D5A", "#310C5A", "#3A0B59", "#4B0A5B", "#51095A",
                "#5B0B5A", "#730861", "#6B015A", "#7C0362", "#820060", "#8B0A62", "#940B63", "#A40563",
                "#AE0964", "#C20F68", "#BD0062", "#C40063", "#C50059", "#C6004B", "#C70542", "#C60A39",
                "#C2224A", "#C50A2A", "#C50B21", "#C5121B", "#C4121A", "#C4232B", "#C3332B", "#C2452A"
        };

        for (size_t i = 0; i < number_colors; ++i) {
            RGBcolors.push_back("\"" + cols[i] + "\"");
        }

        RGBcoeff = (number_colors - 1) / (get_count_genomes() - 1);
    }
}

}
#endif // CONFIG_STRUCT_IMPL_HPP