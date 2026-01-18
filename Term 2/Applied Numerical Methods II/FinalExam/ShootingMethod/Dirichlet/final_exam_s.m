function dy =final_exam_s(x,y)
dy = [y(2); -y(2).^2 - y(1) + log(x)];
end