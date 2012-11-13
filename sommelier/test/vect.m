load in.mat
expected = a * v;
comparison = isalmost(r, expected);
if (!isequal(comparison, ones(size(v))))
  a
  v
  r
  expected
  comparison
  error("Vector multiplication failed.")
endif
