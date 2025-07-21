// Copyright (c) 2025, Davide Cenzato.  All rights reserved.
// Use of this source code is governed
// by a MIT license that can be found in the LICENSE file.

/*
 *  thr_qgram_profile_static_dict: To fill
 */

#ifndef thr_qgram_profile_HPP_
#define thr_qgram_profile_HPP_

#include <tsl/robin_map.h>
#include <ef_compressed_dictionary.hpp>
#include <vector>
#include <limits>

namespace tqd{

template<typename value_type = uint8_t,
		 class static_dictionary = tqd::ef_compressed_dictionary<>>
class thr_qgram_profile{

private:

	static_dictionary profile; // static compressed dictionary mapping q-gram to its occurrence counter
	tsl::robin_map<size_t, value_type> map; // dynamic uncompressed dictionary

	const int q;    // q-gram length
	const int t;    // threshold value
	const int base; // alphabet base

	std::vector<size_t> pows;

public:

	thr_qgram_profile(int q_, int t_, int base_) : q(q_), t(t_), base(base_)
	{
		if(bitsize(t) > sizeof(value_type)*8)
		{
			std::cerr << "The current value type cannot represent the occurrence counter " <<
			             "sizes given threshold = " << t << std::endl;
			exit(1);
		}

		if(exceeds_size_t(base,q))
		{
			std::cerr << "The size of the fingerprints cannot be represented using 64 bits " <<
			             "given base = " << base << " and q-gram length = " << q << std::endl;
			exit(1);
		}

		pows.resize(base);
		size_t pow = base;
		for(size_t i=0;i<(q-2);++i){ pow *= base; }
		for(size_t i=0;i<base;++i){ pows[i] = i*pow; }

		profile.set((pow * base),bitsize(t+1));
	}

	void parse_ASCII_DNAonly_sequence(const std::string& sequence)
	{
		size_t i = 0;
		size_t fingerprint = 0; 

		while((i + q) < sequence.size()) 
		{
			bool init = false;
			while( not init and (i + q) < sequence.size() ) 
			{
				init = init_window_ASCII_DNAonly(fingerprint,i,sequence);
				if( init ){
					add_fingerprint(fingerprint);
				}
				else{ i++; }
			}

			if(init)
				while(i < sequence.size())
				{
					if( not is_DNA[sequence.at(i)] ){ i++; break; }

					fingerprint -= pows[dna_to_code[sequence.at(i-q)]]; // AAAC
					fingerprint *= base;
					fingerprint += dna_to_code[sequence.at(i++)];

					add_fingerprint(fingerprint);
				}
		}
	}

	void compute_compressed_profile()
	{
		if(map.size() > 0)
		{
			std::vector<std::pair<size_t,value_type>> occ_counters;
			occ_counters.resize(map.size());
			size_t i = 0;
		    for(const auto& key_value : map)
		    {
		    	occ_counters[i].first = key_value.first;
		    	occ_counters[i++].second = key_value.second;
		    }
		    map.clear();
		    // sort keys
		    std::sort(occ_counters.begin(), occ_counters.end(), [](auto &left, auto &right)
		    {
			    return left.first < right.first;
			});

			profile.template build<std::vector<std::pair<size_t,value_type>>>(occ_counters);
		}
	}
	
	size_t compute_tqd(const thr_qgram_profile& other, int64_t t_ = -1)
	{
		if(t_ == -1) t_ = this->t;

		size_t n1 = profile.no_keys();
		size_t n2 = other.profile.no_keys();
		size_t distance = 0;

		for(size_t i=0;i<n1;++i)
		{
			auto kv = profile.get_ith_key_value(i+1);

			int64_t val = other.profile.get_value(kv.first);

			if(val == -1){ distance++; }
			else
			{
				val = std::min(val,t_+1);
				int64_t k_val = std::min(static_cast<int64_t>(kv.second),t_+1);

				n2--;
				if(k_val != val){ distance++; }
			}
		}

		return distance + n2;
	}

	size_t serialize(std::ostream& out)
	{
		size_t w_bytes = 0;

		out.write((char*)&q, sizeof(q));
		out.write((char*)&t, sizeof(t));
		out.write((char*)&base, sizeof(base));
		w_bytes += sizeof(q) + sizeof(t) + sizeof(base);

		w_bytes += profile.serialize(out);

		return w_bytes;
	}

	void load(std::istream& in)
	{
		in.read((char*)&q, sizeof(q));
		in.read((char*)&t, sizeof(t));
		in.read((char*)&base, sizeof(base));

		profile.load(in);
	}

private:

	inline void add_fingerprint(size_t fingerprint)
	{
		auto it = this->map.find(fingerprint);
		if(it != this->map.end()){ if(it.value() <= t) it.value() += 1; }
		else{ this->map[fingerprint] = 1; }
	}

	inline bool 
	init_window_ASCII_DNAonly(size_t& fingerprint, size_t& i, const std::string& sequence) const
	{
		fingerprint = 0;
		for(size_t j=0;j<q;++j)
		{
			if( not is_DNA[sequence.at(i)] ){ return false; }
			fingerprint *= base;
			fingerprint += dna_to_code[sequence.at(i++)];
		}

		return true;
	}

	uint8_t bitsize(size_t x)
	{
	    if(x == 0) return 1;
	    return 64 - __builtin_clzll(x);
	}


	bool exceeds_size_t(uint64_t base, uint64_t exp)
	{
	    uint64_t val = 1;
	    for (uint64_t i = 0; i < exp; ++i)
	    {
	        if (val > std::numeric_limits<size_t>::max() / base)
	            return true;
	        val *= base;
	    }
	    return false;
	}

	static constexpr bool is_DNA[128] = {
	    0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,
	    0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,
	    0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,
	    0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,
	    0, 1, 0, 1,  0, 0, 0, 1,  0, 0, 0, 0,  0, 0, 0, 0,
	    0, 0, 0, 0,  1, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,
	    0, 1, 0, 1,  0, 0, 0, 1,  0, 0, 0, 0,  0, 0, 0, 0,
	    0, 0, 0, 0,  1, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0
	};
	static constexpr unsigned char dna_to_code[128] = {
	    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	    4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
	    4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	    4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
	    4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 5
	};

};

}

#endif