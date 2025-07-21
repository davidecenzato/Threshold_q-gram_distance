// Copyright (c) 2025, Davide Cenzato.  All rights reserved.
// Use of this source code is governed
// by a MIT license that can be found in the LICENSE file.

#ifndef UTILS_HPP_
#define UTILS_HPP_

class simple_progress_bar
{
public:

    simple_progress_bar(size_t total, std::string header, size_t width = 70)
        : total(total), width(width), current(0), header(header) {}

    void advance()
    {
        ++current;
        float ratio = current / static_cast<float>(total);
        int64_t c = ratio * width;

        std::cout << "\r\t" << header << " [";
        for (size_t x = 0; x < c; ++x) std::cout << "=";
        if(c < width){ std::cout << ">"; }
        for (size_t x = c+1; x < width; ++x) std::cout << " ";
        std::cout << "] " << static_cast<int64_t>(ratio * 100.0) << "%";
        std::cout.flush();
    }

    void done() { std::cout << "\n" << std::endl; }

private:

    size_t total, width, current;
    std::string header;
};

void create_folder_if_not_exists(const std::string& path)
{
    if (!std::filesystem::exists(path)) {
        std::filesystem::create_directories(path); 
    }
}

std::vector<std::string> format_fasta_string(const std::string& fasta_string)
{
    std::istringstream iss(fasta_string);
    std::string line;
    std::vector<std::string> sequences;
    std::string current_sequence;

    while (std::getline(iss, line))
    {
        if (line.empty()) continue;

        if (line[0] == '>') 
        {
            // new header line
            if (!current_sequence.empty()) {
                sequences.push_back(current_sequence);
                current_sequence.clear();
            }
        } else {
            // append sequence line
            current_sequence += line;
        }
    }

    // add the final sequence
    if (!current_sequence.empty()) {
        sequences.push_back(current_sequence);
    }

    return sequences;
}

std::vector<int> format_int_list(const std::string& input)
{
    std::vector<int> res;
    std::stringstream ss(input);
    std::string token;

    while (std::getline(ss, token, ','))
        { res.push_back(std::stoi(token)); }

    return res;
}

#endif