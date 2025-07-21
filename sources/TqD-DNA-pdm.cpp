#include <iostream>
#include <string>
#include <filesystem>

#include <thr_qgram_profile.hpp>
#include "utils.hpp"

void help(){

    std::cout << "tqd-dna [options]" << std::endl <<
    "Options:" << std::endl <<
    "-h          Print usage info." << std::endl <<
    "-i <arg>    Directory path containing the input FASTA files. (REQUIRED)" << std::endl <<
    "-q <arg>    q-gram lengths list. (Def. 8)" << std::endl <<
    "-t <arg>    Threshold values list. (Def. 1)" << std::endl <<
    "-o <arg>    Output directory path for pairwise distance matrices. (REQUIRED)" << std::endl;
    exit(0);
} 

int main(int argc, char* argv[])
{
    if(argc < 2) { help(); }

    std::string inputPath, outputPath;
    const int base = 4;
    std::string q_list, t_list;

    int opt;
    while ((opt = getopt(argc, argv, "hi:q:t:o:")) != -1)
    {
        switch (opt){
            case 'h':
                help();
            break;
            case 'i':
                inputPath = std::string(optarg);
            break;
            case 'o':
                outputPath = std::string(optarg);
            break;
            case 'q':
            	q_list = std::string(optarg);
            break;
            case 't':
                t_list = std::string(optarg);
            break;
            default:
                help();
            return -1;
        }
    }

    if(inputPath.empty() or outputPath.empty()){ help(); }

    std::vector<int> q_ = format_int_list(q_list);
    std::vector<int> t_ = format_int_list(t_list);
    int max_t = 1;
    for(const int& t : t_){ max_t = std::max(max_t,t); }

    std::cout << "\n[INFO] Computing the pairwise distance matrices of the sequences" 
              << " in " << inputPath << " for the following parameters: q = [" << q_list 
              << "] t = [" << t_list << "]" << std::endl;

    std::cout << "\n[INFO] Computing the Threshold q-gram profiles for the files" 
              << " in " << inputPath << std::endl << std::endl;

    size_t no_files = 0;
    for (const auto& filepath : std::filesystem::directory_iterator(inputPath)){ no_files++; }

    create_folder_if_not_exists(outputPath);

    for(const auto& q : q_)
    { // TqD profiles computation
        int t = max_t;
        std::string qfolder = outputPath+"/q_"+std::to_string(q);
        create_folder_if_not_exists(qfolder);

        std::string header = "\t[q = " + std::to_string(q) + ", t = " + std::to_string(t) +  "] ";
        simple_progress_bar bar(no_files,header);

        for (const auto& filepath : std::filesystem::directory_iterator(inputPath))
        {
            std::string file = filepath.path();
            std::string input_sequence;

            std::ifstream in(file,std::ifstream::binary);
            in.seekg(0, std::ifstream::end);
            input_sequence.resize(in.tellg());
            in.seekg(0, std::ifstream::beg);

            in.read(static_cast<char *>(&input_sequence[0]),input_sequence.size());
            in.close();
            auto cleaned_sequences = format_fasta_string(input_sequence);

            tqd::thr_qgram_profile<> profile(q,t,base);
            for(const auto& s : cleaned_sequences)
                profile.parse_ASCII_DNAonly_sequence(s);
            
            profile.compute_compressed_profile();

            std::string output_file = qfolder + "/" + std::string(filepath.path().filename()) + ".tqd";
            std::ofstream o(output_file.c_str(),std::ofstream::binary);
            profile.serialize(o);
            o.close();

            bar.advance();
        }

        bar.done();
    }

    std::cout << "\n[INFO] Computing the Threshold q-gram distances for the profiles" 
              << " of sequences in " << inputPath << " and storing the pairwise distance matrices"
              << "in " << outputPath << std::endl << std::endl;

    for(const auto& q : q_)
        for(const auto& t : t_)
        {
            std::string qfolder = outputPath+"/q_"+std::to_string(q);

            size_t no_files = 0;
            std::vector<std::string> file_list;
            std::vector<std::string> name_list;
            for (const auto& filepath : std::filesystem::directory_iterator(qfolder))
            { 
                no_files++; 
                file_list.push_back(filepath.path());
                std::string name = filepath.path().filename();
                name_list.push_back(name.substr(0,name.find('.')));
            }

            auto sum_up_to_n = [](int n) { return n * (n - 1) / 2; };
            size_t no_combinations = sum_up_to_n(no_files);

            std::string header = "\t[q = " + std::to_string(q) + ", t = " + std::to_string(t) +  "] ";
            simple_progress_bar bar(no_combinations,header);

            std::string outFile    = outputPath+"/q_"+std::to_string(q)+"t_"+std::to_string(t)+".tsv";
            std::ofstream out(outFile,std::ifstream::binary);

            std::vector<uint64_t> distance_list;
            distance_list.reserve(no_combinations);

            for(size_t i=0;i<file_list.size()-1;++i)
            {
                std::ifstream in(file_list[i],std::ifstream::binary);
                tqd::thr_qgram_profile<> profile1(q,t,base);
                profile1.load(in);
                in.close();

                for(size_t j=i+1;j<file_list.size();++j)
                {
                    std::ifstream in_(file_list[j],std::ifstream::binary);
                    tqd::thr_qgram_profile<> profile2(q,t,base);
                    profile2.load(in_);
                    in_.close();

                    auto distance =  profile1.compute_tqd(profile2,t);
                    out << name_list[i] << "\t" << name_list[j] << "\t" << distance << "\n";
                    distance_list.push_back(distance);

                    bar.advance();
                }
            }

            bar.done();
            out.close();
        }

    std::cout << "\n[INFO] Deleting temporary profile files" << std::endl;

    for(const auto& q : q_)
    {
        std::string qfolder = outputPath+"/q_"+std::to_string(q);
        std::filesystem::remove_all(qfolder);
    }

    std::cout << "\n[INFO] Done!" << std::endl;

    return 0;
}