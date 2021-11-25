// Utility

// LIST

function list_to_array(L) {
  const A = [];
  let i = 0;
  for (let p = L; !is_null(p); p = tail(p)) {
    A[i] = head(p);
    i = i + 1;
  }
  return A;
}

function multi_map(f, xss) {
  if (is_null(head(xss))) {
    return null;
  } else {
    return pair(f(map(head, xss)), multi_map(f, map(tail, xss)));
  }
}

// Destructive functions
function d_filter(pred, xs) {
  if (is_null(xs)) {
    return xs;
  } else if (pred(head(xs))) {
    set_tail(xs, d_filter(pred, tail(xs)));
    return xs;
  } else {
    return d_filter(pred, tail(xs));
  }
}

// Check list of nums all different
function all_different(nums) {
  if (is_null(nums)) {
    return true;
  } else {
    const head_is_unique = is_null(member(head(nums), tail(nums)));
    return head_is_unique && all_different(tail(nums));
  }
}

// Pascal / Combinations
function pascal(row, pos) {
  return pos === 1 || pos === row
    ? 1
    : pascal(row - 1, pos - 1) + pascal(row - 1, pos);
}

// ARRAY
function swap(A, i, j) {
  const temp = A[i];
  A[i] = A[j];
  A[j] = temp;
}

function copy_array(A) {
  const len = array_length(A);
  const B = [];
  for (let i = 0; i < len; i = i + 1) {
    B[i] = A[i];
  }
  return B;
}

function reverse_array(A) {
  const len = array_length(A);
  const half_len = math_floor(len / 2);
  for (let i = 0; i < half_len; i = i + 1) {
    swap(A, i, len - 1 - i);
  }
}

function array_to_list(A) {
  const len = array_length(A);
  let L = null;
  for (let i = len - 1; i >= 0; i = i - 1) {
    L = pair(A[i], L);
  }
  return L;
}

function map_array(f, arr) {
  const len = array_length(arr);
  function iter(i) {
    if (i < len) {
      arr[i] = f(arr[i]);
      iter(i + 1);
    }
  }
  iter(0);
}

function accumulate_array(op, init, A) {
  let ans = copy_array(init);
  const length = array_length(A);

  function helper(i) {
    if (i < length) {
      ans = op(ans, A[i]);
      helper(i + 1);
    }
  }
  helper(0);

  return ans;
}

function filter_array(pred, A) {
  // YOUR SOLUTION HERE
  const length = array_length(A);
  const ans = [];

  function helper(i) {
    if (i < length) {
      if (pred(A[i])) {
        ans[array_length(ans)] = A[i];
      }
      helper(i + 1);
    }
  }
  helper(0);

  return ans;
}

// STREAMS
function add_streams(s1, s2) {
  return is_null(s1)
    ? s2
    : is_null(s2)
    ? s1
    : pair(head(s1) + head(s2), () =>
        add_streams(stream_tail(s1), stream_tail(s2))
      );
}

function scale_stream(c, stream) {
  return stream_map((x) => c * x, stream);
}

function stream_pairs(s) {
  return is_null(s)
    ? null
    : stream_append(
        stream_map((sn) => pair(head(s), sn), stream_tail(s)),
        stream_pairs(stream_tail(s))
      );
}

function mul_stream(s1, s2) {
  return pair(head(s1) * head(s2), () =>
    add_streams(
      stream_tail(scale_stream(head(s2), s1)),
      mul_stream(stream_tail(s2), s1)
    )
  );
}
// TREES
function scale_tree(tree, factor) {
  return map(
    (sub_tree) =>
      !is_list(sub_tree) ? factor * sub_tree : scale_tree(sub_tree, factor),
    tree
  );
}

function map_tree(f, tree) {
  return map(
    (sub_tree) => (!is_list(sub_tree) ? f(sub_tree) : map_tree(f, sub_tree)),
    tree
  );
}

function count_data_items(tree) {
  return is_null(tree)
    ? 0
    : (is_list(head(tree)) ? count_data_items(head(tree)) : 1) +
        count_data_items(tail(tree));
}

// BAE
function evaluate_BAE_tree(bae_tree) {
  if (is_list(bae_tree)) {
    const left = evaluate_BAE_tree(head(bae_tree));
    const right = evaluate_BAE_tree(head(tail(tail(bae_tree))));
    const op = head(tail(bae_tree));
    if (op === "+") {
      return left + right;
    } else if (op === "-") {
      return left - right;
    } else if (op === "*") {
      return left * right;
    } else {
      // (op === "/")
      return left / right;
    }
  } else {
    // is a number
    return bae_tree;
  }
}

function build_BAE_tree(bae_list) {
  let next_token = bae_list;

  function build_tree() {
    if (equal(head(next_token), "(")) {
      next_token = tail(next_token);
      const left_tree = build_tree();
      const op = head(next_token);
      next_token = tail(next_token);
      const right_tree = build_tree();
      next_token = tail(next_token); // skip over ")"
      return list(left_tree, op, right_tree);
    } else {
      // token is a number
      const token = head(next_token);
      next_token = tail(next_token);
      return token;
    }
  }

  return build_tree();
}

// BST
// A binary search tree (BST) is either null or a list with three elements: a number x, a left BST,
// and a right BST, where every number in the left BST is smaller than the number x, and every
// number in the right BST is larger than the number x.
function BST_min(bst) {
  return is_null(bst)
    ? Infinity
    : is_null(head(tail(bst)))
    ? head(bst)
    : BST_min(head(tail(bst)));
}

function BST_find(x, bst) {
  return is_null(bst)
    ? false
    : x === head(bst)
    ? true
    : x < head(bst)
    ? BST_find(x, head(tail(bst)))
    : BST_find(x, head(tail(tail(bst))));
}

function BST_to_list(bst) {
  if (is_null(bst)) {
    return null;
  } else {
    const ltree = head(tail(bst));
    const num = head(bst);
    const rtree = head(tail(tail(bst)));
    return append(BST_to_list(ltree), pair(num, BST_to_list(rtree)));
  }
}

function accumulate_bst(op, initial, bst) {
  // SOLUTION 1:
  if (is_empty_binary_tree(bst)) {
    return initial;
  } else {
    const s = accumulate_bst(op, initial, right_subtree_of(bst));
    const t = op(value_of(bst), s);
    return accumulate_bst(op, t, left_subtree_of(bst));
  }
  // // SOLUTION 2:
  // function listify_bst(b) {
  //   if (is_empty_binary_tree(b)) {
  //     return null;
  //   } else {
  //     const left_list = listify_bst(left_subtree_of(b));
  //     const value = value_of(b);
  //     const right_list = listify_bst(right_subtree_of(b));
  //     return append(left_list, pair(value, right_list));
  // }
  // }
  // return accumulate(op, initial, listify_bst(bst));
}

// MEMOIZATION
function read(n, k) {
  return is_undefined(mem[n]) ? undefined : mem[n][k];
}

function write(n, k, value) {
  if (is_undefined(mem[n])) {
    mem[n] = [];
  }
  mem[n][k] = value;
}

function memoize(f) {
  const mem = [];
  function mf(x) {
    if (mem[x] !== undefined) {
      return mem[x];
    } else {
      const result = f(x);
      mem[x] = result;
      return result;
    }
  }
  return mf;
}

function mchoose(n, k) {
  if (read(n, k) !== undefined) {
    return read(n, k);
  } else {
    const result =
      k > n
        ? 0
        : k === 0 || k === n
        ? 1
        : mchoose(n - 1, k) + mchoose(n - 1, k - 1);
    write(n, k, result);
    return result;
  }
}

// ACTIVE LISTS
function act_append(as, bs) {
  const len_as = act_length(as);
  return (pos) => (pos < len_as ? as(pos) : bs(pos - len_as));
}

// SEARCHING
function binary_search_array(A, v) {
  let low = 0;
  let high = array_length(A) - 1;
  while (low <= high) {
    const mid = math_floor((low + high) / 2);
    if (v === A[mid]) {
      break;
    } else if (v < A[mid]) {
      high = mid - 1;
    } else {
      low = mid + 1;
    }
  }
  return low <= high;
}

// SORTING
function bubblesort_array(A) {
  const len = array_length(A);
  for (let i = len - 1; i >= 1; i = i - 1) {
    for (let j = 0; j < i; j = j + 1) {
      if (A[j] > A[j + 1]) {
        const temp = A[j];
        A[j] = A[j + 1];
        A[j + 1] = temp;
      }
    }
  }
}

function bubblesort_list(L) {
  const len = length(L);
  for (let i = len - 1; i >= 1; i = i - 1) {
    let p = L;
    for (let j = 0; j < i; j = j + 1) {
      if (head(p) > head(tail(p))) {
        const temp = head(p);
        set_head(p, head(tail(p)));
        set_head(tail(p), temp);
      }
      p = tail(p);
    }
  }
}

function selection_sort_list(xs) {
  if (is_null(xs)) {
    return xs;
  } else {
    const x = smallest(xs);
    return pair(x, selection_sort_list(remove(x, xs)));
  }
}

function selection_sort_array(A) {
  const len = array_length(A);
  for (let i = 0; i < len - 1; i = i + 1) {
    let min_pos = find_min_pos(A, i, len - 1);
    swap(A, i, min_pos);
  }
}
function find_min_pos(A, low, high) {
  let min_pos = low;
  for (let j = low + 1; j <= high; j = j + 1) {
    if (A[j] < A[min_pos]) {
      min_pos = j;
    }
  }
  return min_pos;
}

function insertion_sort(A) {
  const len = array_length(A);
  for (let i = 1; i < len; i = i + 1) {
    let j = i - 1;
    while (j >= 0 && A[j] > A[j + 1]) {
      swap(A, j, j + 1);
      j = j - 1;
    }
  }
}

function take(xs, n) {
  return n === 0 ? null : pair(head(xs), take(tail(xs), n - 1));
}

function drop(xs, n) {
  return n === 0 ? xs : drop(tail(xs), n - 1);
}

function merge_sort(xs) {
  if (is_null(xs) || is_null(tail(xs))) {
    return xs;
  } else {
    const mid = middle(length(xs));
    return merge(merge_sort(take(xs, mid)), merge_sort(drop(xs, mid)));
  }
}

function merge_sort_takedrop(xs) {
  if (is_null(xs) || is_null(tail(xs))) {
    return xs;
  } else {
    const td = take_drop(xs, middle(length(xs)));
    return merge(merge_sort_takedrop(head(td)), merge_sort_takedrop(tail(td)));
  }
}

function merge_sort_array(A) {
  merge_sort_helper(A, 0, array_length(A) - 1);
}

function merge_sort_helper(A, low, high) {
  if (low < high) {
    const mid = math_floor((low + high) / 2);
    merge_sort_helper(A, low, mid);
    merge_sort_helper(A, mid + 1, high);
    merge(A, low, mid, high);
  }
}

function merge(A, low, mid, high) {
  const B = []; // temporary array
  let left = low;
  let right = mid + 1;
  let Bidx = 0;
  while (left <= mid && right <= high) {
    if (A[left] <= A[right]) {
      B[Bidx] = A[left];
      left = left + 1;
    } else {
      B[Bidx] = A[right];
      right = right + 1;
    }
    Bidx = Bidx + 1;
  }
  while (left <= mid) {
    B[Bidx] = A[left];
    Bidx = Bidx + 1;
    left = left + 1;
  }
  while (right <= high) {
    B[Bidx] = A[right];
    Bidx = Bidx + 1;
    right = right + 1;
  }
  for (let k = 0; k < high - low + 1; k = k + 1) {
    A[low + k] = B[k];
  }
}

function take_drop(xs, n) {
  function helper(ys, k, acc) {
    return k === 0
      ? pair(acc, ys)
      : helper(tail(ys), k - 1, pair(head(ys), acc));
  }
  return helper(xs, n, null);
}

function partition(xs, p) {
  // Your answer from Task 1
  return pair(
    filter((a) => a <= p, xs),
    filter((a) => a > p, xs)
  );
}

function quicksort(xs) {
  // Your answer here
  return is_null(xs)
    ? null
    : is_null(tail(xs))
    ? xs
    : accumulate(
        append,
        null,
        list(
          quicksort(head(partition(tail(xs), head(xs)))),
          list(head(xs)),
          quicksort(tail(partition(tail(xs), head(xs))))
        )
      );
}

// MISC
function get_sublist(start, end, L) {
  function helper(pos, ys, result) {
    // SOLUTION 1:
    if (pos < start) {
      return helper(pos + 1, tail(ys), result);
    } else if (pos <= end) {
      return helper(pos + 1, tail(ys), pair(head(ys), result));
    } else {
      return reverse(result);
    }
    // SOLUTION 2:
    // return pos < start
    //   ? helper(pos + 1, tail(ys), result)
    //   : pos <= end
    //   ? helper(pos + 1, tail(ys), pair(head(ys), result))
    //   : reverse(result);
  }
  return helper(0, L, null);
}

function insertions(x, ys) {
  return map(
    (k) => append(take(ys, k), pair(x, drop(ys, k))),
    enum_list(0, length(ys))
  );
}

function permutations(ys) {
  return is_null(ys)
    ? list(null)
    : accumulate(
        append,
        null,
        map((x) => map((p) => pair(x, p), permutations(remove(x, ys))), ys)
      );
}

function are_permutation(xs1, xs2) {
  return is_null(xs1) && is_null(xs2)
    ? true
    : !is_null(xs1) && !is_null(xs2)
    ? !is_null(member(head(xs1), xs2)) &&
      are_permutation(tail(xs1), remove(head(xs1), xs2))
    : false;
}

function remove_duplicates(lst) {
  return accumulate(
    (x, xs) => (is_null(member(x, xs)) ? pair(x, xs) : xs),
    null,
    lst
  );
}

// Generate all possible subsets
function subsets(xs) {
  return accumulate(
    (x, ss) =>
      append(
        ss,
        map((s) => pair(x, s), ss)
      ),
    list(null),
    xs
  );
}

function is_subset(S, T) {
  if (is_null(S)) {
    return true;
  } else if (is_null(T)) {
    return false;
  } else if (head(S) < head(T)) {
    return false;
  } else if (head(S) === head(T)) {
    return is_subset(tail(S), tail(T));
  } else {
    return is_subset(S, tail(T));
  }
}

function combinations(xs, k) {
  if (k === 0) {
    return list(null);
  } else if (is_null(xs)) {
    return null;
  } else {
    const s1 = combinations(tail(xs), k - 1);
    const s2 = combinations(tail(xs), k);
    const x = head(xs);
    const has_x = map((s) => pair(x, s), s1);
    return append(has_x, s2);
  }
}

function check_parentheses(paren_list) {
  function check(count, xs) {
    if (is_null(xs)) {
      return count === 0;
    } else if (count < 0) {
      return false;
    } else if (head(xs) === "(") {
      return check(count + 1, tail(xs));
    } else {
      // (head(xs) === ")")
      return check(count - 1, tail(xs));
    }
  }

  return check(0, paren_list);
}

function sum_of_list(xs) {
  function sum(still_to_process, sum_so_far) {
    if (is_null(still_to_process)) {
      return sum_so_far;
    } else {
      return sum(tail(still_to_process), sum_so_far + head(still_to_process));
    }
  }
  return sum_of_list(xs, 0);
}

function rotate_matrix(M) {
  const n = array_length(M); // M is assumed n x n.
  function swap(r1, c1, r2, c2) {
    const temp = M[r1][c1];
    M[r1][c1] = M[r2][c2];
    M[r2][c2] = temp;
  }
  // Do a matrix transpose first.
  for (let r = 0; r < n; r = r + 1) {
    for (let c = r + 1; c < n; c = c + 1) {
      swap(r, c, c, r);
    }
  }
  // Then reverse each row.
  const half_n = math_floor(n / 2);
  for (let r = 0; r < n; r = r + 1) {
    for (let c = 0; c < half_n; c = c + 1) {
      swap(r, c, r, n - c - 1);
    }
  }
}

function count_pairs(x) {
  let pairs = null;
  function check(y) {
    if (!is_pair(y)) {
      return undefined;
    } else if (!is_null(member(y, pairs))) {
      return undefined;
    } else {
      pairs = pair(y, pairs);
      check(head(y));
      check(tail(y));
    }
  }
  check(x);
  return length(pairs);
}

function tetrate(b, n) {
  return n === 1 ? b : math_pow(b, tetrate(b, n - 1));
}

function fib_list(N) {
  function helper(count, rev) {
    return count === N
      ? rev
      : helper(count + 1, pair(head(rev) + head(tail(rev)), rev));
  }
  return reverse(helper(2, list(1, 0)));
}

// Sorts the array of numbers in ascending order.
function sort_ascending(A) {
  const len = array_length(A);
  for (let i = 1; i < len; i = i + 1) {
    const x = A[i];
    let j = i - 1;
    while (j >= 0 && A[j] > x) {
      A[j + 1] = A[j];
      j = j - 1;
    }
    A[j + 1] = x;
  }
}

function digits_to_string(digits) {
  const len = array_length(digits);
  let str = "";
  for (let i = 0; i < len; i = i + 1) {
    str = str + stringify(digits[i]);
  }
  return str;
}

function distinct(xs) {
  if (is_null(xs) || is_null(tail(xs))) {
    return true;
  } else {
    if (is_null(member(head(xs), tail(xs)))) {
      return distinct(tail(xs));
    } else {
      return false;
    }
  }
}
