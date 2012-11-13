load in.mat
expected = a * alpha;
comparison = isalmost(r, expected);
if (!isequal(comparison, ones(size(a))))
  a
  alpha
  r
  expected
  comparison
  error("Array scaling failed.")
endif
