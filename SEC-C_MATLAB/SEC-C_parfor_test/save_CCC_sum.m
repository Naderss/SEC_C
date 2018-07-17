%save CCC_sum in CCC_dir directory 

function save_CCC_sum(CCC_sum,CCC_dir,i)
save([CCC_dir,'CCC_sum_',num2str(i),'.mat'],'CCC_sum');
end

