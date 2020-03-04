function [param, index] = getParamIndex(i)
  if (i<=16*5)
    j = i;
    param = ceil(j/16);
    index = mod(j-1, 16)+1;
  elseif (i<=16*5+32*5)
    j = i-16*5;
    param = ceil(j/32) + 5;
    index = mod(j-1, 32)+1;
  elseif (i<=16*5+32*5+48*5)
    j = i-16*5-32*5;
    param = ceil(j/48) + 10;
    index = mod(j-1, 48)+1;
  elseif (i<=16*5+32*5+48*5+64*5)
    j = i-16*5-32*5-48*5;
    param = ceil(j/64) + 15;
    index = mod(j-1, 64)+1;
  end
end