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
	template <typename T, typename Key>
	class Tree
	{
	protected:
		Node<T>* root;
		size_t count;
		virtual bool _find(Node<T>* subroot, const Key& key, T& item) const = 0;
		virtual Node<T>* _insert(Node<T>* subroot, const T& item) = 0;
		virtual Node<T>* _remove(Node<T>* subroot , const Key& key, Node<T>*& item) = 0;
		virtual void _print(Node<T>* subroot, int level) const = 0;
		virtual Node<T>* removeMin(Node<T>* subroot, Node<T>*& min) = 0;
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
		bool remove(const Key& key, T& item) {
			Node<T>* t = NULL;
			root = _remove(root, key, t);
			// key is not found
			if (t == NULL) return false;
			item = t->value();
			--count;
			delete t;
			return true;
		}
		bool removeAny(T& item) {
			// delete minimum value
			if (root == NULL) return false;
			Node<T>* t;
			root = removeMin(root, t);
			item = t->value();
			delete t;
			--count;
			return true;
		}
		bool find(const Key& key, T& item) const {
			return _find(root, key, item);
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
	class BSTree : public Tree<T, T>
	{
	private:
		bool _find(Node<T>* subroot, const T& key, T& item) const {
			if (subroot == NULL) {
				return false;
			} else if (key < subroot->value) {
				// check left
				return _find(subroot->left, key, item);
			} else if (key > subroot->value) {
				// check right
				return _find(subroot->left, key, item);
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
		// return the subtree after the node with the specified key has been removed
		// remove designated node by replacing it with the node
		//  with the least key value greater than the one being removed,
		//  i.e. the node with the least key value in the right branch
		//  (alternatively, can replace the node with the greatest key value
		//  less than the one being removed, i.e. in the left branch)
		Node<T>* _remove(Node<T>* subroot , const T& key, Node<T>*& item) {
			if (subroot == NULL) {
				// item is not in tree
				return NULL;
			} else if (key < subroot->value) {
				// check left
				subroot->left = _remove(subroot->left, key, item);
			} else if (key > subroot->value) {
				subroot->right = _remove(subroot->right, key, item);
			} else {
				// found it: remove it
				Node<T>* tmp;
				item = subroot;
				if (subroot->left == NULL) {
					// only a right child: point to it
					subroot = subroot->right;
				} else if (subroot->right == NULL) {
					// only a left child: point to it
					subroot = subroot->left;
				} else {
					// node to be removed has both children
					// find node with least key value that is greater than the node;
					//  store it in $tmp
					// (no need to to assignment?)
					subroot->right = removeMin(subroot->right, tmp);
					// swap values
					T value = subroot->value;
					subroot->value = tmp->value;
					tmp->value = value;
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
		
		// return the subroot with the min key value node removed
		//  also returns pointer to the removed node by reference
		Node<T>* removeMin(Node<T>* subroot, Node<T>*& min) {
			if (subroot->left == NULL) {
				// node with min key has no left-child
				// (if it did, it would not have the min key)
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
	
	// T must support T.keys()
	// the size of T.keys() must match ndim
	template <typename T>
	class KDTree : public Tree<T, int*>
	{
		typedef int* Keys;
	private:
		size_t ndim;
		bool _find(Node<T>* subroot, const Keys& coord, T& item) const {
			return this->_find(subroot, coord, item, 0);
		}
		bool _find(Node<T>* subroot, const Keys& coord, T& item, int discrim) const {
			if (subroot == NULL) return false;  // empty tree
			Keys currCoord = subroot->value.keys();
			bool coordsAreEqual = true;
			for (size_t i = 0; i < ndim; ++i) {
				if (currCoord[i] != coord[i]) {
					coordsAreEqual = false;
					break;
				}
			}
			if (coordsAreEqual) {
				item = subroot->value;
				return true;
			}
			if (coord[discrim] < currCoord[discrim]) {
				// search left branch, incrementing the discriminant dimension
				return _find(subroot->left, coord, item, (discrim+1)%ndim);
			} else {
				// search righ branch
				return _find(subroot->right, coord, item, (discrim+1)%ndim);
			}
		}
		Node<T>* _insert(Node<T>* subroot, const T& item) {
			return this->_insert(subroot, item, 0);
		}
		Node<T>* _insert(Node<T>* subroot, const T& item, int discrim) {
			if (subroot == NULL) {
				// empty tree: create node
				return new Node<T>(item);
			}
			int* currCoord = subroot->value.keys();
			int* coord = item.keys();
			if (coord[discrim] < currCoord[discrim]) {
				// insert on left
				subroot->left = _insert(subroot->left, item, (discrim+1)%ndim);
			} else {
				// insert on right
				subroot->right = _insert(subroot->right, item, (discrim+1)%ndim);
			}
			// return subtree with node inserted
			return subroot;
		}
		void _print(Node<T>* subroot, int level) const {
			if (subroot == NULL) return;  // empty tree: nothing to print
			_print(subroot->left, level+1);
			// indent to level
			for (int i = 0; i < level; ++i) {
				cout << "- ";
			}
			// print tuple
			cout << "/";
			for (int i = 0; i < ndim; ++i) {
				cout << subroot->value.keys()[i] << "/";
			}
			cout << endl;
			_print(subroot->right, level+1);
		}
		// Remving node N
		// first find node N, then replace its record by the record in N's right
		//   subtree with the least value of N's discriminator, or
		//   by the record in N's left subtree with the greatest value for this
		//   discriminator
		Node<T>* _remove(Node<T>* subroot, const Keys& key, Node<T>*& item) {
			throw runtime_error("KDTree::_remove(...) is yet not implemented!");
		}
		Node<T>* removeMin(Node<T>* subroot, Node<T>*& min) {
			throw runtime_error("KDTree::_removeMin(...) is yet not implemented!");
		}
		// auxilary function for removing a node
		// finds the node with the mininum key value for the specified discriminator
		Node<T>* findMin(Node<T>* subroot, int discrim, int currDiscrim) {
			if (subroot == NULL) return NULL;
			Node<T> *a, *b;
			Keys *keys, *akeys, *bkeys;
			keys = subroot->value.keys();
			a = findMin(subroot->left, discrim, (currDiscrim+1)%ndim);
			if (a != NULL) akeys = a->value.keys();
			if (discrim != currDiscrim) {
				// discriminator for the current level is different from the
				//   the specified discrminator
				// Thus, min could be on either side: check the right side too
				// Otherwise, having checked the left side sufficient
				b = findMin(subroot->right, discrim, (currDiscrim+1)%ndim);
				if (b != NULL) bkeys = b->value.keys();
				if (a == NULL || (b != NULL && bkeys[discrim] < akeys[discrim])) {
					// Right side has a smaller key value
					a = b;
					akeys = bkeys;
				}
			}
			// Now, a has the smallest value in children
			if (a == NULL || (keys[discrim] < akeys[discrim])) {
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