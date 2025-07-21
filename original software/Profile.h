// structure for a node of the RB-tree
typedef struct structure{
        unsigned long long int qgram; // number that represent the q-gram
        int size1; // number of this qgram in the first strings
        int size2; // number of this qgram in the second strings
        struct structure* sx;
        struct structure* dx;
        struct structure* father;
        bool color; //0 is black and 1 is red
        }node;


// class that represent the threshold-qgram profile
// We represent it using an RB-tree
class Profile{
	public:
		//variables
		int threshold;
		
		//constructor
		Profile(int thr);
		//destructor
		~Profile();
		
		//insert a node in the tree
		void insert(unsigned long long int qgram, int n_str);
		
		//calculate the distance
		unsigned long long int calculateDistance();
		
		//get the total status
		int* get_total_status();
		
		//get the number of nodes
		unsigned long long int getNumber_of_nodes();
		
	private:
		node *root;
		node *nil; //special node that represent NULL and is black
		bool update; //said if the total_status is update, every time we calculate the distance is update
		int *total_status;
		unsigned long long int number_of_nodes;
		
		unsigned long long int calculateDistanceRecursive(node *n);
		
		//support for insert
		void left_rotate(node *x);
		void right_rotate(node *x);
		void RB_insert_fixup(node *z);
};
