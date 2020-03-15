function saveArray(fn,title,A);
	io = open(fn,"w");
	for name in title
		write(io,name);
		write(io,";");
	end
	write(io,"\n");
	for i in 1:size(A)[1]
		for k in 1:size(A)[2]
			write(io,string(A[i,k]));
			write(io,";");
		end
		write(io,"\n");
	end
	close(io);
end
