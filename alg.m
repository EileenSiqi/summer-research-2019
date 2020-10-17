%microbial genetic algorithm
P=20; %population
D=5; %deme
rec=0.9; %recombination rate
mut=0.9; %mutation rate
mut_id=0.04; %individual mutation rate
vrange=0.4;
epsilonrange=5.2;
alpharange=3.2;
tour=12000;
num=6; %num of gene
gene=zeros(P,num);
for i=1:P
    v=-0.2+vrange*rand(2,1);
    epsilon=epsilonrange*rand(2,1);
    alpha=-1.6+alpharange*rand(2,1); %generate a gene
    gene(i,:)=cat(1,v,epsilon,alpha); 
    f(i)=fitness(gene(i,:));
end
for i=1:tour
    A=ceil(P*rand());
    B=mod((A+ceil(D*rand())),P)+1;
    if f(A)<f(B)
        W=A; L=B;
    else
        W=B; L=A;
    end
    if rand()<rec %do the recombination process
        for j=1:num
            if rand()<2/num
                gene(L,j)=gene(W,j);
            end
        end
    end
    if rand()<mut %do the mutation process
        mut_pos=ceil(num*rand());
        if mut_pos/num<=1/3
            gene(L,mut_pos)=gene(L,mut_pos)+vrange*mut_id*randn();
        elseif mut_pos/num<=2/3
            gene(L,mut_pos)=gene(L,mut_pos)+epsilonrange*mut_id*randn();
        else
            gene(L,mut_pos)=gene(L,mut_pos)+alpharange*mut_id*randn();
        end
    end
    f(L)=fitness(gene(L,:));
    Min(i)=min(f);
    Mean(i)=mean(f);
end
plot([1:tour],Mean,[1:tour],Min);
xlabel('tournament number');
ylabel('fitness');
