function fieldArray = getFieldArray(S, fieldName, rowNum)
% S 2-D struct Array [exploreParamNo, exploreTimes]
% fieldName only support 1-D vector or scalar fields
arrayLength = length(S(rowNum,:));
if (~isscalar(S(rowNum,1).(fieldName)))
  fieldArray = zeros(arrayLength, length(S(rowNum,1).(fieldName)));
  for i=1:arrayLength
    fieldArray(i,:) = S(rowNum, i).(fieldName);
  end
else
  fieldArray = zeros(1, arrayLength);
  for i=1:arrayLength
    fieldArray(i) = S(rowNum, i).(fieldName);
  end
end
