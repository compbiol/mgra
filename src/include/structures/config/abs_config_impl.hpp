//
// Created by pavel on 1/8/16.
//

#ifndef MGRA_ABS_CONFIG_IMPL_HPP
#define MGRA_ABS_CONFIG_IMPL_HPP

namespace config {

/**
 * Load functions
 */
inline void abs_config::load_genomes(Json::Value const &genomes) {
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

inline void abs_config::load_genome(Json::Value const &genome, size_t index) {
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
            if (!alias.isString()) {
                std::cerr << "ERROR: Alias must have string type" << std::endl;
                exit(1);
            }

            if (genome_number.count(alias.asString()) > 0) {
                std::cerr << "ERROR: Duplicate alias " << alias.asString() << std::endl;
                exit(1);
            }
            genome_number.insert(std::make_pair(alias.asString(), index));
        }
    }
}

inline void abs_config::load_block_type(Json::Value const &type_of_file) {
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

inline void abs_config::load_files(Json::Value const &path_to_files) {
    if (!path_to_files.isArray()) {
        std::cerr << "ERROR: files section is array with at least one blocks file" << std::endl;
        exit(1);
    }

    std::string dir_with_cfg = path::parent_path(config_file_path);

    for (auto const & path : path_to_files) {
        if (!path.isString()) {
            std::cerr << "ERROR: path to file must have string type" << std::endl;
            exit(1);
        }
        path_to_blocks_file.push_back(path::make_full_path(path::append_path(dir_with_cfg, path.asString())));
    }
}

inline void abs_config::load_tree(Json::Value const &tree) {
    if (!tree.isString()) {
        std::cerr << "ERROR: tree must have string type" << std::endl;
        exit(1);
    }

    using namespace structure::phyl_tree;
    phylotree = TreeBuilder<phylogeny_tree_t>::build_tree(tree.asString());
}

inline void abs_config::load_output_directory(Json::Value const &path_to_dir) {
    if (!path_to_dir.isString()) {
        std::cerr << "ERROR: path_to_dir field must have string type" << std::endl;
        exit(1);
    }

    out_path_directory = path_to_dir.asString();
}
/**
 * Save functions
 */
inline Json::Value abs_config::save_genomes() const {
    Json::Value genomes(Json::arrayValue);
    for (size_t ind = 0; ind != priority_name.size(); ++ind) {
        genomes.append(save_genome(ind));
    }
    return genomes;
}

inline Json::Value abs_config::save_genome(size_t ind) const {
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

inline Json::Value abs_config::save_block_type() const {
    Json::Value type;

    if (block_file_type == block_file_type_t::infercars) {
        type = Json::Value("infercars");
    } else if (block_file_type == block_file_type_t::grimm) {
        type = Json::Value("grimm");
    }

    return type;
}

inline Json::Value abs_config::save_files() const {
    Json::Value files(Json::arrayValue);

    for (auto const & path : path_to_blocks_file) {
        files.append(Json::Value(path));
    }

    return files;
}

inline Json::Value abs_config::save_tree() const {
    std::ostringstream out;
    structure::phyl_tree::NewickTreePrinter<phylogeny_tree_t> printer(out);
    printer.print_tree(phylotree);
    Json::Value tree(out.str());
    return tree;
}

inline Json::Value abs_config::save_output_directory() const {
    return Json::Value(out_path_directory);
}

}

#endif //MGRA_ABS_CONFIG_IMPL_HPP
