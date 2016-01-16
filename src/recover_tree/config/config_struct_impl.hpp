//
// Created by pavel on 10/21/15.
//

#ifndef MGRA_CONFIG_STRUCT_IMPL_HPP
#define MGRA_CONFIG_STRUCT_IMPL_HPP

/*
 * Read configure file in JSON format and init all config
 */
template<class mcolor_t>
void load(tree_config<mcolor_t> &cfg, std::istream & source) {
    INFO("Start load cfg file in MGRA format");

    Json::Value root;
    Json::Reader reader;
    bool isError =  reader.parse(source, root);

    if (!isError) {
        std::string error_msg = "Failed to parse JSON config\n" + reader.getFormattedErrorMessages();
        std::cerr << error_msg << std::endl;
        exit(1);
    }

    cfg.load(root);
    INFO("Finish load cfg file in MGRA format")
}

/*
 * Different function, which load our input.
 */
template<class mcolor_t>
void tree_config<mcolor_t>::load(Json::Value const &root) {
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

    if (!root["trees"].isNull()) {
        load_trees(root["trees"]);
    }

    if (!root["statistics"].isNull()) {
        load_source_statistics(root["statistics"]);
    } else {
        recover_tree_statistic = distribution;
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

    default_rgb_colors();
}

template <class mcolor_t>
Json::Value tree_config<mcolor_t>::save() const {
    Json::Value root;

    root["genomes"] = save_genomes();

    root["format_of_blocks"] = save_block_type();

    root["files"] = save_files();

    root["trees"] = save_trees();

    root["output_directory"] = save_output_directory();

    root["statistics"] = save_source_statistics();

    return root;
}

/**
 * Load functions
 */
template<class mcolor_t>
void tree_config<mcolor_t>::load_genomes(Json::Value const &genomes) {
    if (!genomes.isArray()) {
        std::cerr << "ERROR: genomes section is array with at least three genomes" << std::endl;
        exit(1);
    }

    priority_name.resize(genomes.size());
    size_t i = 0;

    for (auto const & genome : genomes) {
        load_genome(genome, i++);
    }
}

template<class mcolor_t>
void tree_config<mcolor_t>::load_genome(Json::Value const &genome, size_t index) {
    if (!genome.isObject()) {
        std::cerr << "ERROR: genome section must be object" << std::endl;
        exit(1);
    }

    // Parse optional section about genome id
    if (!genome["genome_id"].isNull()) {
        if (!genome["genome_id"].isInt()) {
            std::cerr << "ERROR: genome_id field in genome section must have int type" << std::endl;
            exit(1);
        }
        index = (size_t) genome["genome_id"].asInt();
    }

    // Parse required section priority name
    {
        if (genome["priority_name"].isNull()) {
            std::cerr << "ERROR: priority name is required field for genome section" << std::endl;
            exit(1);
        }

        if (!genome["priority_name"].isString()) {
            std::cerr << "ERROR: priority_name field in genome section must have string type" << std::endl;
            exit(1);
        }

        auto const &name = genome["priority_name"].asString();
        if (genome_number.count(name) > 0) {
            std::cerr << "ERROR: Genome identificator " << name << " is not unique!" << std::endl;
            exit(1);
        }

        priority_name[index] = name;
        genome_number.insert(std::make_pair(name, index));
        number_to_genome.insert(std::make_pair(index, name));
    }

    // Parse optional section about genome aliases names
    if (!genome["aliases"].isNull()) {
        if (!genome["aliases"].isArray()) {
            std::cerr << "ERROR: aliases field in genomes section is array with at least one names" << std::endl;
            exit(1);
        }
        for (auto const & alias : genome["aliases"]) {
            if (genome_number.count(alias.asString()) > 0) {
                std::cerr << "ERROR: Duplicate alias " << alias.asString() << std::endl;
                exit(1);
            }
            genome_number.insert(std::make_pair(alias.asString(), index));
        }
    }
}

template<class mcolor_t>
void tree_config<mcolor_t>::load_block_type(Json::Value const &type_of_file) {
    if (!type_of_file.isString()) {
        std::cerr << "ERROR: type format of block files must have string format" << std::endl;
        exit(1);
    }

    if (type_of_file.asString() == "infercars") {
        block_file_type = block_file_type_t::infercars;
    } else if (type_of_file.asString() == "grimm") {
        block_file_type = block_file_type_t::grimm;
    } else {
        std::cerr << "ERROR: Unsuported type of blocks (grimm, infercars required)" << std::endl;
    }
}

template<class mcolor_t>
void tree_config<mcolor_t>::load_files(Json::Value const &path_to_files) {
    if (!path_to_files.isArray()) {
        std::cerr << "ERROR: files section is array with at least one blocks file" << std::endl;
        exit(1);
    }

    std::string dir_with_cfg = path::parent_path(config_file_path);

    for (auto const & file : path_to_files) {
        if (!file.isString()) {
            std::cerr << "ERROR: path to file must have string type" << std::endl;
            exit(1);
        }
        path_to_blocks_file.push_back(path::make_full_path(path::append_path(dir_with_cfg, file.asString())));
    }
}

template<class mcolor_t>
void tree_config<mcolor_t>::load_trees(Json::Value const &trees) {
    if (!trees.isArray()) {
        std::cerr << "ERROR: trees section is array with one full tree or at least one subtrees" << std::endl;
        exit(1);
    }

    for (auto const & tree : trees) {
        load_tree(tree);
    }
}

template<class mcolor_t>
void tree_config<mcolor_t>::load_tree(Json::Value const &tree) {
    if (!tree.isString()) {
        std::cerr << "ERROR: tree must have string type" << std::endl;
        exit(1);
    }

    using namespace structure::phyl_tree;
    phylotrees.push_back(TreeBuilder<phylogeny_tree_t>::build_tree(tree.asString()));
}

template<class mcolor_t>
void tree_config<mcolor_t>::load_output_directory(Json::Value const &path_to_dir) {
    if (!path_to_dir.isString()) {
        std::cerr << "ERROR: path_to_dir field must have string type" << std::endl;
        exit(1);
    }

    out_path_directory = path_to_dir.asString();
}

template <class mcolor_t>
void tree_config<mcolor_t>::load_source_statistics(Json::Value const &source_stat) {
    if (!source_stat.isString()) {
        std::cerr << "ERROR: source_statistic field must have string type" << std::endl;
        exit(1);
    }

    if (source_stat.asString() == "patterns") {
        recover_tree_statistic = patterns;
    } else if (source_stat.asString() == "distribution") {
        recover_tree_statistic = distribution;
    }
}

/**
 * Save functions
 */
template<class mcolor_t>
Json::Value tree_config<mcolor_t>::save_genomes() const {
    Json::Value genomes(Json::arrayValue);
    for (size_t ind = 0; ind != priority_name.size(); ++ind) {
        genomes.append(save_genome(ind));
    }
    return genomes;
}

template<class mcolor_t>
Json::Value tree_config<mcolor_t>::save_genome(size_t ind) const {
    Json::Value genome(Json::objectValue);

    genome["genome_id"] = Json::Value((int) ind);

    genome["priorety_name"] = Json::Value(priority_name[ind]);

    genome["aliases"] = Json::Value(Json::arrayValue);
    auto range = number_to_genome.equal_range(ind);
    for (auto it = range.first; it != range.second; ++it) {
        genome["aliases"].append(Json::Value(it->second));
    }

    return genome;
}

template<class mcolor_t>
Json::Value tree_config<mcolor_t>::save_block_type() const {
    Json::Value type;

    if (block_file_type == block_file_type_t::infercars) {
        type = Json::Value("infercars");
    } else if (block_file_type == block_file_type_t::grimm) {
        type = Json::Value("grimm");
    }

    return type;
}

template<class mcolor_t>
Json::Value tree_config<mcolor_t>::save_files() const {
    Json::Value files(Json::arrayValue);

    for (auto const & file : path_to_blocks_file) {
        files.append(Json::Value(file));
    }

    return files;
}

template<class mcolor_t>
Json::Value tree_config<mcolor_t>::save_trees() const {
    Json::Value trees(Json::arrayValue);

    for (size_t ind = 0; ind != phylotrees.size(); ++ind) {
        trees.append(save_tree(ind));
    }

    return trees;
}

template<class mcolor_t>
Json::Value tree_config<mcolor_t>::save_tree(size_t ind) const {
    assert(ind != phylotrees.size());
    std::ostringstream out;
    //writer::TXT_NewickTree<phylogeny_tree_t> printer(out);
    //printer.print_tree(phylotrees[ind]); //FIXME
    Json::Value tree(out.str());
    return tree;
}

template<class mcolor_t>
Json::Value tree_config<mcolor_t>::save_output_directory() const {
    return Json::Value(out_path_directory);
}

template<class mcolor_t>
Json::Value tree_config<mcolor_t>::save_source_statistics() const {
    if (recover_tree_statistic == patterns) {
        return Json::Value("patterns");
    } else {
        return Json::Value("distribution");
    }
}

/**
* Function, which initilization structure on default parameters
*/
template<class mcolor_t>
void tree_config<mcolor_t>::default_directory_organization() {
    std::string const LOGGER_FILENAME = "mgra.log";
    std::string const TREES_FILENAME = "trees";
    std::string const NEWICK_FILENAME = "summary.newick";

    out_path_to_logger_file = path::append_path(out_path_directory, LOGGER_FILENAME);
    trees_path = path::append_path(out_path_directory, TREES_FILENAME);
    tree_summary_path = path::append_path(out_path_directory, NEWICK_FILENAME);
}

template<class mcolor_t>
void tree_config<mcolor_t>::default_rgb_colors() {
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


#endif //MGRA_CONFIG_STRUCT_IMPL_HPP
