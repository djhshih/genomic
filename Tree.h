#ifndef genomic_Tree_h
#define genomic_Tree_h

#include <iostream>
#include <stdexcept>

using namespace std;

namespace tree {
	
	template <typename T>
	struct Node
	{
		T value;
		Node<T> *left, *right;
		bool leaf;
		Node() {}
		Node(T nodeValue)
		: value(nodeValue), leaf(true), left(NULL), right(NULL) {}
		Node(T nodeValue, Node* leftChild, Node* rightChild)
		: value(nodeValue), left(leftChild), right(rightChild), leaf(false) {}
	};
	
	// Tree derived classes have same interface, different mechanism
	// Use Public Function Calls Virtual Protected idiom
	template <typename T>
	class Tree
	{
	protected:
		Node<T>* root;
		size_t count;
		virtual bool _find(Node<T>* subroot, const T& value, T& item) const = 0;
		virtual Node<T>* _insert(Node<T>* subroot, const T& item) = 0;
		virtual Node<T>* _remove(Node<T>* subroot , const T& value, Node<T>*& item) = 0;
		virtual void _print(Node<T>* subroot, int level) const = 0;
	private:
		void __clear(Node<T>* subroot) {
			if (subroot == NULL) return;
			__clear(subroot->left);
			__clear(subroot->right);
			delete subroot;
		}
	public:
		Tree() : root(NULL), count(0) {}
		virtual ~Tree() {
			clear();
		}
		void clear() {
			__clear(root);
			root = NULL;
			count = 0;
		}
		void insert(const T& item) {
			root = _insert(root, item);
			++count;
		}
		bool remove(const T& value, T& item) {
			Node<T>* t = NULL;
			root = _remove(root, value, t);
			// value is not found
			if (t == NULL) return false;
			item = t->value();
			--count;
			delete t;
			return true;
		}
		bool find(const T& value, T& item) const {
			return _find(root, value, item);
		}
		unsigned int size() {
			return count;
		}
		void print() const {
			if (root == NULL) {
				cout << "The tree is empty." << endl;
			} else {
				_print(root, 0);
			}
		}
	};
	
	template <typename T>
	class BSTree : public Tree<T>
	{
	private:
		bool _find(Node<T>* subroot, const T& value, T& item) const {
			if (subroot == NULL) {
				return false;
			} else if (value < subroot->value) {
				// check left
				return _find(subroot->left, value, item);
			} else if (value > subroot->value) {
				// check right
				return _find(subroot->left, value, item);
			} else {
				// found
				item = subroot->value;
				return true;
			}
		}
		Node<T>* _insert(Node<T>* subroot, const T& item) {
			if (subroot == NULL) {
				// empty tree: create node
				return new Node<T>(item);
			}
			if (item < subroot->value) {
				// insert on left
				subroot->left = _insert(subroot->left, item);
			} else {
				// insert on right
				subroot->right = _insert(subroot->right, item);
			}
			// return subtree with node inserted
			return subroot;
		}
		// return the subtree after the node with the specified value has been removed
		// remove designated node by replacing it with the node
		//  with the least value value greater than the one being removed,
		//  i.e. the node with the least value value in the right branch
		//  (alternatively, one can replace the node with the greatest value value
		//  less than the one being removed, i.e. in the left branch; however,
		//  duplicated values would result in an invalid tree)
		Node<T>* _remove(Node<T>* subroot , const T& value, Node<T>*& item) {
			if (subroot == NULL) {
				// item is not in tree
				return NULL;
			} else if (value < subroot->value) {
				// check left
				subroot->left = _remove(subroot->left, value, item);
			} else if (value > subroot->value) {
				subroot->right = _remove(subroot->right, value, item);
			} else {
				// found it: remove it
				Node<T>* tmp;
				item = subroot;
				if (subroot->left == NULL) {
					// only a right child (or no child): point to it
					subroot = subroot->right;
				} else if (subroot->right == NULL) {
					// only a left child: point to it
					subroot = subroot->left;
				} else {
					// Node to be removed has both children
					// Replace current node's value with the node storing the
					//  least value value that is greater than the node;
					//  store it in $tmp (passed by reference)
					subroot->right = removeMin(subroot->right, tmp);
					// swap values
					T value = subroot->value;
					subroot->value = tmp->value;
					tmp->value = value;
					// Now, subroot is the replacing node
					//   and $tmp stores the removed node: pass it by reference
					item = tmp;
				}
			}
			return subroot;
		}
		// print tree with the root at the left, branching to the right
		void _print(Node<T>* subroot, int level) const {
			if (subroot == NULL) return;  // empty tree: nothing to print
			_print(subroot->left, level+1);
			// indent to level
			for (int i = 0; i < level; ++i) {
				cout << "- ";
			}
			cout << subroot->value << endl;
			_print(subroot->right, level+1);
		}
		
		// return the subroot with the min value value node removed
		//  also returns pointer to the removed node by reference
		Node<T>* removeMin(Node<T>* subroot, Node<T>*& min) {
			if (subroot->left == NULL) {
				// node with min value has no left-child
				// (if it did, it would not have the min value)
				min = subroot;
				// past the right child back to the parent,
				//  s.t. the parent can assign the left child to it,
				//  effectively removing $min from the tree
				return subroot->right;
			} else {
				// continue traversing down the left branch to find the min
				subroot->left = removeMin(subroot->left, min);
				return subroot;
			}
		}
	public:
		BSTree() {}
		~BSTree() {}
	};
	
	
	
	template <typename T>
	class Container
	{
	public:
		typedef int Key;
		typedef int* Keys;
	private:
		const size_t ndim;
	public:
		Keys keys;
		T value;
		Container(size_t nKeys) : ndim(nKeys), keys(new Key[nKeys]) {}
		Container(const Container& d)
		: ndim(d.ndim), keys(new Key[d.ndim]), value(d.value) {
			for (size_t i = 0; i < d.ndim; ++i) {
				keys[i] = d.keys[i];
			}
		}
		~Container() {
			delete keys;
		}
		Container& operator=(const Container& c) {
			Container tmp(c);
			swap(this->keys, tmp.keys);
			return *this;
		}
		bool operator==(const Container& c) const {
			return (*this) == c.keys;
		}
		bool operator==(const Keys& k) const {
			bool equal = true;
			for (size_t i = 0; i < ndim; ++i) {
				if (keys[i] != k[i]) {
					equal = false;
					break;
				}
			}
			return equal;
		}
		Key& operator[](size_t i) const {
			return keys[i];
		}
	};
	
	// KDKeyedData is an class with a k-dimensional key set
	//   The number of keys must match KDTree::ndim
	//   (Is there an efficient way to enforce this?)
	// Required methods:
	//   KDKeyedData::operator==(KDKeyedData&) for comparing keys of two KDKeyedData objects
	//   KDKeyedData::operater[](size_t) for accessing a key in the key set
  // Container is a class that fulfills these requirements
	template <typename KDKeyedData>
	class KDTree : public Tree<KDKeyedData>
	{
	private:
		size_t ndim;
		bool _find(Node<KDKeyedData>* subroot, const KDKeyedData& value, KDKeyedData& item) const {
			return this->_find(subroot, value, item, 0);
		}
		bool _find(Node<KDKeyedData>* subroot, const KDKeyedData& value, KDKeyedData& item, int discrim) const {
			if (subroot == NULL) return false;  // empty tree
			if (subroot->value == value) {
				item = subroot->value;
				return true;
			}
			if (value[discrim] < subroot->value[discrim]) {
				// search left branch, incrementing the discriminant dimension
				return _find(subroot->left, value, item, (discrim+1)%ndim);
			} else {
				// search righ branch
				return _find(subroot->right, value, item, (discrim+1)%ndim);
			}
		}
		Node<KDKeyedData>* _insert(Node<KDKeyedData>* subroot, const KDKeyedData& item) {
			return this->_insert(subroot, item, 0);
		}
		Node<KDKeyedData>* _insert(Node<KDKeyedData>* subroot, const KDKeyedData& item, int discrim) {
			if (subroot == NULL) {
				// empty tree: create node
				return new Node<KDKeyedData>(item);
			}
			if (subroot->value[discrim] < item[discrim]) {
				// insert on left
				subroot->left = _insert(subroot->left, item, (discrim+1)%ndim);
			} else {
				// insert on right
				subroot->right = _insert(subroot->right, item, (discrim+1)%ndim);
			}
			// return subtree with node inserted
			return subroot;
		}
		void _print(Node<KDKeyedData>* subroot, int level) const {
			if (subroot == NULL) return;  // empty tree: nothing to print
			_print(subroot->left, level+1);
			// indent to level
			for (int i = 0; i < level; ++i) {
				cout << "- ";
			}
			// print tuple
			cout << "/";
			for (int i = 0; i < ndim; ++i) {
				cout << subroot->value[i] << "/";
			}
			cout << endl;
			_print(subroot->right, level+1);
		}
		
		Node<KDKeyedData>* _remove(Node<KDKeyedData>* subroot, const KDKeyedData& value, Node<KDKeyedData>*& item) {
			return _remove(subroot, value, item, 0);
		}
		// Remving node N
		// If N has one child...
		//   We cannot simply assign N's parent to point to N's child, since
		//   the lefts would be changed, create an invalid KD tree.
		//   Therefore, replace N by an appropriate node in the left or right tree.
		// If N has only a right child (or both children)...
		//   First find node N, then replace its record by the record in N's right
		//   subtree with the least value of N's discriminator
		// If the right subtree does not exist, it is not satisfactory to
		//   replace N's record with the recording having the greatest value for
		//   the discriminator in the left subtree: this value might appear more than
		//   once in the left subtree.
		//   This would lead to the equal values for the discriminator in N's left tree.
		//   (Recall that nodes in the left tree must have lower values of the
		//   discriminator.)
		// Solution: first move the left subtree of node N to become the right
		//   subtree. Then, process with normal deletion process, replacing the
		//   the record of N to be deleted with the record containing the
		//   *least* value of the discriminator from what is now N's right subtree
		Node<KDKeyedData>* _remove(Node<KDKeyedData>* subroot, const KDKeyedData& value, Node<KDKeyedData>*& item, int discrim) {
			
			if (subroot == NULL) {
				// item is not in tree
				return NULL;
			}
			
			if (subroot->value == value) {
				// Found N == subroot: proceed to remove it
				
				if (subroot->right == NULL && subroot->left == NULL) {
					// N has no children: N can simply be removed
					// Pass N by reference back to caller
					item = subroot;
					// Set subroot to NULL, s.t. calling function can effectively 
					//   remove N from tree while setting its pointer to subroot
					subroot = NULL;
				} else {
					// N has at least one child
				
					// Always replace N's record by the record in N's right subtree with
					//   the least value of N's discriminator
					
					// One special case: right subtree does not exist
					//   Move to the left subtree to the right subtree
					if (subroot->right == NULL) {
						// only a left child: point to it
						subroot->right = subroot->left;
						subroot->left = NULL;
					}
					// Find the node with the minimum key of the N's discriminator,
					//   in the right subtree
					//   $discrim == N's discrminator, since subroot == N
					Node<KDKeyedData>* min = findMin(subroot->right, discrim, discrim);
					// At this point, N must have a right child, because if it only
					//   had a left child, it has been assigned as the right child
					// Since N has a right child, min cannot be NULL
					// After finding min, remove it from the right subtree
					subroot->right = _remove(subroot->right, min->value, min, discrim);
					// swap values between N and the min Node
					KDKeyedData value = subroot->value;
					subroot->value = min->value;
					min->value = value;
					// Now, min holds of the record of N
					item = min;
				}
				
			} else if (value[discrim] < subroot->value[discrim]) {
				// search left branch, incrementing the discriminant dimension
				// assign left child to the residual left tree with the target node removed
				subroot->left = _remove(subroot->left, value, item, (discrim+1)%ndim);
			} else {
				// search righ branch
				// assign right child to the residual right tree with the target node removed
				subroot->right = _remove(subroot->right, value, item, (discrim+1)%ndim);
			}
			
			return subroot;
		}
		// Finds and returns the node with the mininum key value of the
		//   specified discriminator
		// Auxilary function for removing a node
		Node<KDKeyedData>* findMin(Node<KDKeyedData>* subroot, int discrim, int currDiscrim) {
			if (subroot == NULL) return NULL;
			Node<KDKeyedData> *a, *b;
			a = findMin(subroot->left, discrim, (currDiscrim+1)%ndim);
			if (discrim != currDiscrim) {
				// discriminator for the current level is different from the
				//   the specified discrminator
				// Thus, min could be on either side: check the right side too
				// Otherwise, having checked the left side sufficient
				b = findMin(subroot->right, discrim, (currDiscrim+1)%ndim);
				if (a == NULL || (b != NULL && b->value[discrim] < a->value[discrim])) {
					// Right side has a smaller key value
					a = b;
				}
			}
			// Now, a has the smallest value in children
			if (a == NULL || (subroot->value[discrim] < a->value[discrim])) {
				return subroot;
			}
			return a;
		}
	public:
		KDTree(size_t nDimensions) : ndim(nDimensions) {}
		~KDTree() {}
	};
	
};


#endif