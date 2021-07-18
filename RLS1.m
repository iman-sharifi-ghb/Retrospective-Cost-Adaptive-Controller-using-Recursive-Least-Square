function  [teta , P] = RLS1(phi , y  , teta , P , Nc)
K = P * phi * (eye(size(y)) + phi' * P * phi)^-1 ;
teta = teta + K * (y - phi' * teta) ;
P = (eye(Nc) - K * phi')*P ;
end