// Copyright (c) 2025, Davide Cenzato.  All rights reserved.
// Use of this source code is governed
// by a MIT license that can be found in the LICENSE file.

#ifndef ef_compressed_dictionary_HPP_
#define ef_compressed_dictionary_HPP_

#include <sdsl/int_vector.hpp>
#include <sdsl/sd_vector.hpp>

#include <chrono>

namespace tqd{

template<class compressed_bitvector = sdsl::sd_vector<>, 
		 class select_key_support   = sdsl::select_support_sd<>,
		 class rank_key_support     = sdsl::rank_support_sd<>,
	     class bitpacked_vector     = sdsl::int_vector<>>
class ef_compressed_dictionary{

private:

	compressed_bitvector keys; // compressed bitvector containing the keys
	select_key_support select; // select keys support
	rank_key_support rank; // rank keys support
	bitpacked_vector values;   // bitpacked vector containing the values

	size_t u; // universe size
	size_t n; // number of keys
	uint8_t w; // values width

public:

	ef_compressed_dictionary(){}

	void set(size_t u_, uint8_t w_)
	{
		u = u_;
		w = w_;
		values.width(w);
	}

	template<typename occ_type>
	void build(const occ_type occ_counters)
	{
		sdsl::sd_vector_builder builder(u,occ_counters.size());
		for(auto idx: occ_counters){ builder.set(idx.first); }
		keys = sdsl::sd_vector<>(builder);
		// set bitvector len
		n = occ_counters.size();

		values.resize(n);
		size_t i = 0;
		for(auto idx: occ_counters){ values[i++] = idx.second; }
	}
	
	size_t no_keys() const { return n; }

	std::pair<size_t,size_t> get_ith_key_value(size_t i) const
	{
		assert(i < n+1);

		size_t key = select(i);

		return std::make_pair(key,values[rank(key+1)-1]);
	}

	int64_t get_value(size_t key) const
	{
		assert(key < u);

		if(keys[key])
		{
			size_t r = rank(key+1);
			return values[r-1];
		}
		else{ return -1; }
	}

	size_t serialize(std::ostream& out)
	{
		size_t w_bytes = 0;

		out.write((char*)&u, sizeof(u));
		out.write((char*)&n, sizeof(n));
		out.write((char*)&w, sizeof(w));
		w_bytes += sizeof(u) + sizeof(n) + sizeof(w);

		w_bytes += keys.serialize(out);
		sdsl::util::init_support(select,&keys);
		sdsl::util::init_support(rank,&keys);
		w_bytes += values.serialize(out);

		return w_bytes;
	}

	void load(std::istream& in)
	{
		in.read((char*)&u, sizeof(u));
		in.read((char*)&n, sizeof(n));
		in.read((char*)&w, sizeof(w));

		keys.load(in);
		sdsl::util::init_support(select,&keys);
		sdsl::util::init_support(rank,&keys);
		values.load(in);
	}

};

}

#endif