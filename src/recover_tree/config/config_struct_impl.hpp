//
// Created by pavel on 10/21/15.
//

#ifndef MGRA_CONFIG_STRUCT_IMPL_HPP
#define MGRA_CONFIG_STRUCT_IMPL_HPP

template <class mcolor_t>
void load(tree_config<mcolor_t>& cfg, std::string const& filename) {
    INFO("Start load cfg file for recover tree")

    Json::Value root;
    Json::Reader reader;
    std::ifstream input(filename);
    bool isError =  reader.parse(input, root);

    if (!isError) {
        std::string error_msg = "Failed to parse JSON config\n" + reader.getFormattedErrorMessages();
        std::cerr << error_msg << std::endl;
        exit(1);
    }

    cfg.load(root);

    INFO("Finish load cfg file for recover tree")
}

template <class mcolor_t>
void tree_config<mcolor_t>::load(Json::Value const & root) {
    if (!root["genomes"].isNull()) {
        load_genomes(root["genomes"]);
    } else {
        std::cerr << "ERROR: genomes section is required" << std::endl;
        exit(1);
    }

    if (!root["files"].isNull()) {
        load_files(root["files"]);
    } else {
        std::cerr << "ERROR: files section is required" << std::endl;
        exit(1);
    }

    if (!root["output_directories"].isNull()) {
        ;
    }

    if (!root["trees"].isNull()) {
        load_trees(root["trees"]);
    }

}

template <class mcolor_t>
void tree_config<mcolor_t>::load_genomes(Json::Value const& genomes) {
    if (!genomes.isArray()) {
        std::cerr << "ERROR: genomes section is array with at least three genomes" << std::endl;
        exit(1);
    }

    priority_name.resize(genomes.size());
    size_t i = 0;
    for (Json::Value::iterator it = genomes.begin(); it != genomes.end(); ++it) {
        load_genome(*it, i++);
    }
}

template <class mcolor_t>
void tree_config<mcolor_t>::load_genome(Json::Value const& genome, size_t index) {
    if (!genome["genome_id"].isNull()) {

        if (!genome["genome_id"].isInt()) {
            std::cerr << "ERROR: genome_id field in genome section must have int type" << std::endl;
            exit(1);
        }

        index = genome["genome_id"].asInt();
    }

    if (!genome["priority_name"].isNull()) {

        if (!genome["priority_name"].isString()) {
            std::cerr << "ERROR: priority_name field in genome section must have string type" << std::endl;
            exit(1);
        }

        auto const & name = genome["priority_name"].asString();
        if (genome_number.count(name) > 0) {
            std::cerr << "ERROR: Genome identificator " << name << " is not unique!" << std::endl;
            exit(1);
        }

        priority_name[index] = name;
        genome_number.insert(std::make_pair(name, index));
    } else {
        std::cerr << "ERROR: priority name is required field for genome section" << std::endl;
        exit(1);
    }

    if (!genome["aliases"].isNull()) {

        if (!genome["aliases"].isArray()) {
            std::cerr << "ERROR: aliases field in genomes section is array with at least one names" << std::endl;
            exit(1);
        }

        for (Json::Value::iterator alias = genome["aliases"].begin(); alias != genome["aliases"].end(); ++alias) {
            if (genome_number.count(alias->asString()) > 0) {
                std::cerr << "ERROR: Duplicate alias " << alias->asString() << std::endl;
                exit(1);
            }
            genome_number.insert(std::make_pair(alias->asString(), index));
        }
    }
}

template <class mcolor_t>
void tree_config<mcolor_t>::load_files(Json::Value const& path_to_files) {
    if (!path_to_files.isArray()) {
        std::cerr << "ERROR: files section is array with at least one blocks file" << std::endl;
        exit(1);
    }

    for (Json::Value::iterator it = path_to_files.begin(); it != path_to_files.end(); ++it) {
        if (!it->isString()) {
            std::cerr << "ERROR: path to file must have string type" << std::endl;
            exit(1);
        }
        paths_to_blocks_file.push_back(it->asString());
    }
}

template <class mcolor_t>
void tree_config<mcolor_t>::load_trees(Json::Value const& trees) {
    if (!trees.isArray()) {
        std::cerr << "ERROR: trees section is array with one full tree or at least one subtrees" << std::endl;
        exit(1);
    }

    for (Json::Value::iterator it = trees.begin(); it != trees.end(); ++it) {
        load_tree(*it);
    }
}

template <class mcolor_t>
void tree_config<mcolor_t>::load_tree(Json::Value const& tree) {
    if (!tree.isString()) {
        std::cerr << "ERROR: tree must have string type" << std::endl;
        exit(1);
    }

    using tree_t = typename structure::BinaryTree<mcolor_t>;
    using node_t = typename tree_t::colored_node_t;

    std::function<void(std::shared_ptr<const node_t>)> get_names_lambda = [&](std::shared_ptr<const node_t> current) -> void {
        if (current->get_left_child()) {
            get_names_lambda(current->get_left_child());
        }

        mcolor_name.insert(std::make_pair(current->get_data(), current->get_name()));

        if (current->get_right_child()) {
            get_names_lambda(current->get_right_child());
        }
    };

    phylotrees.push_back(phylogeny_tree_t(tree.asString(), genome_number, priority_name));
    get_names_lambda(phylotrees.crbegin()->get_root());
}


#endif //MGRA_CONFIG_STRUCT_IMPL_HPP
