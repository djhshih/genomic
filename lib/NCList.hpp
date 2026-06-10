#ifndef cna_NCList_h
#define cna_NCList_h

#include <algorithm>
#include <iterator>
#include <stdexcept>
#include <vector>

#include "typedefs.h"

namespace cna {

class NCList {
public:
	struct IntervalRef {
		position start;
		position end;
		size_t id;
	};

	NCList()
	: size_(0)
	{}

	explicit NCList(const std::vector<IntervalRef>& intervals)
	: size_(0)
	{
		build(intervals);
	}

	void clear()
	{
		starts_.clear();
		ends_.clear();
		source_ids_.clear();
		root_.children.clear();
		root_.ids.clear();
		size_ = 0;
	}

	bool empty() const
	{
		return size_ == 0;
	}

	size_t size() const
	{
		return size_;
	}

	void build(const std::vector<IntervalRef>& intervals)
	{
		clear();
		size_ = intervals.size();
		starts_.reserve(size_);
		ends_.reserve(size_);
		source_ids_.reserve(size_);

		std::vector<size_t> order;
		order.reserve(size_);
		for (size_t i = 0; i < intervals.size(); ++i) {
			if (intervals[i].start > intervals[i].end)
				throw std::logic_error("NCList::build(): start > end");
			starts_.push_back(intervals[i].start);
			ends_.push_back(intervals[i].end);
			source_ids_.push_back(intervals[i].id);
			order.push_back(i);
		}

		std::sort(order.begin(), order.end(), Comparator(*this));

		std::vector<BuildFrame> stack;
		stack.reserve(size_);
		for (size_t oi = 0; oi < order.size(); ++oi) {
			size_t id = order[oi];
			while (!stack.empty() && ends_[stack.back().internal_id] < ends_[id])
				stack.pop_back();

			Node* landing = stack.empty() ? &root_ : stack.back().node;
			landing->ids.push_back(id);
			landing->children.push_back(Node());
			Node* child = &landing->children.back();
			stack.push_back(BuildFrame(child, id));
		}
	}

	template <typename OutputIt>
	void findOverlaps(position start, position end, OutputIt out) const
	{
		forEachOverlap(start, end, OutputReporter<OutputIt>(out));
	}

	bool overlapsAny(position start, position end) const
	{
		return forEachOverlap(start, end, AnyReporter());
	}

	struct InvariantSummary {
		bool sibling_starts_monotone;
		bool sibling_ends_monotone;
		bool parent_contains_children;
	};

	InvariantSummary checkInvariants() const
	{
		InvariantSummary summary;
		summary.sibling_starts_monotone = true;
		summary.sibling_ends_monotone = true;
		summary.parent_contains_children = true;
		checkNodeInvariants(root_, npos, summary);
		return summary;
	}

private:
	struct Node {
		std::vector<Node> children;
		std::vector<size_t> ids;
	};

	struct BuildFrame {
		Node* node;
		size_t internal_id;
		BuildFrame(Node* n, size_t id) : node(n), internal_id(id) {}
	};

	struct WalkFrame {
		const Node* parent;
		size_t child_index;
		WalkFrame(const Node* p, size_t i) : parent(p), child_index(i) {}
	};

	template <typename OutputIt>
	struct OutputReporter {
		mutable OutputIt out;
		OutputReporter(OutputIt o) : out(o) {}
		bool operator()(size_t id) const
		{
			*out++ = id;
			return false;
		}
	};

	struct AnyReporter {
		bool operator()(size_t) const
		{
			return true;
		}
	};

	struct Comparator {
		const NCList& self;
		Comparator(const NCList& s) : self(s) {}
		bool operator()(size_t a, size_t b) const
		{
			if (self.starts_[a] != self.starts_[b])
				return self.starts_[a] < self.starts_[b];
			if (self.ends_[a] != self.ends_[b])
				return self.ends_[a] > self.ends_[b];
			return a < b;
		}
	};

	static const size_t npos = static_cast<size_t>(-1);

	std::vector<position> starts_;
	std::vector<position> ends_;
	std::vector<size_t> source_ids_;
	Node root_;
	size_t size_;

	size_t findLandingChild(const Node& node, position min_end) const
	{
		if (node.ids.empty())
			return npos;
		size_t left = 0;
		size_t right = node.ids.size();
		while (left < right) {
			size_t mid = left + (right - left) / 2;
			if (ends_[node.ids[mid]] < min_end) {
				left = mid + 1;
			} else {
				right = mid;
			}
		}
		return left == node.ids.size() ? npos : left;
	}

	const Node* moveToChild(std::vector<WalkFrame>& stack,
				       const Node* parent,
				       size_t child_index) const
	{
		stack.push_back(WalkFrame(parent, child_index));
		return &parent->children[child_index];
	}

	const Node* moveToRightSiblingOrUncle(std::vector<WalkFrame>& stack) const
	{
		while (!stack.empty()) {
			WalkFrame& top = stack.back();
			if (++top.child_index < top.parent->ids.size())
				return &top.parent->children[top.child_index];
			stack.pop_back();
		}
		return NULL;
	}

	const Node* moveToRightUncle(std::vector<WalkFrame>& stack) const
	{
		if (stack.empty())
			return NULL;
		stack.pop_back();
		return moveToRightSiblingOrUncle(stack);
	}

	template <typename Reporter>
	bool forEachOverlap(position start, position end, Reporter reporter) const
	{
		if (start > end)
			throw std::logic_error("NCList overlap query: start > end");
		if (empty())
			return false;

		const size_t n = findLandingChild(root_, start);
		if (n == npos)
			return false;

		std::vector<WalkFrame> stack;
		stack.reserve(16);
		const Node* node = moveToChild(stack, &root_, n);
		while (node != NULL) {
			WalkFrame& frame = stack.back();
			size_t rgid = frame.parent->ids[frame.child_index];
			if (starts_[rgid] > end) {
				node = moveToRightUncle(stack);
				continue;
			}
			if (reporter(source_ids_[rgid]))
				return true;
			size_t child_index = findLandingChild(*node, start);
			node = child_index != npos ?
				moveToChild(stack, node, child_index) :
				moveToRightSiblingOrUncle(stack);
		}
		return false;
	}

	void checkNodeInvariants(const Node& node, size_t parent_id, InvariantSummary& summary) const
	{
		for (size_t i = 0; i < node.ids.size(); ++i) {
			size_t id = node.ids[i];
			if (i > 0) {
				size_t prev = node.ids[i - 1];
				if (starts_[prev] > starts_[id])
					summary.sibling_starts_monotone = false;
				if (ends_[prev] > ends_[id])
					summary.sibling_ends_monotone = false;
			}
			if (parent_id != npos) {
				if (starts_[id] < starts_[parent_id] || ends_[id] > ends_[parent_id])
					summary.parent_contains_children = false;
			}
			checkNodeInvariants(node.children[i], id, summary);
		}
	}
};

} // namespace cna

#endif
