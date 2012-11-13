load in.mat
expected = a * b;
comparison = isalmost(r, expected);
if (!isequal(comparison, ones(size(a))))
  a
  b
  r
  expected
  comparison
  error("Matrix multiplication failed.")
endif
