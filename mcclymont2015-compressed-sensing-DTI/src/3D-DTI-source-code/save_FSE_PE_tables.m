close all
clear all
clc

mkdir REDACTED
cd REDACTED

etl = 8;


for nv = 64
    
    
    num_slices = nv;
    num_PE = nv;
    
    for accel_factor = 1

        disp(accel_factor)

        % reset the random stream so that we can get the same tables again if needed
        s = RandStream('mcg16807','Seed',0);
        RandStream.setGlobalStream(s);

        for rep = 1:30


            central_region = 0.15;

            [mask, PE1, PE2] = generate_FSE_DTI_mask(accel_factor, ...
                central_region, num_slices, num_PE, etl);

            file_name = ['PE_' num2str(num_PE) 'x' num2str(num_slices) '_etl' num2str(etl) '_f' num2str(accel_factor)  '_rep' num2str(rep) '_DCFill' num2str(central_region) '_pe1'];

            % write the first file
            fid = fopen(file_name, 'w');
            str = 't1 =';
            fprintf(fid,'%s\n',str);
            for i = 1:size(PE1,1)
                seq = PE1(i,:);
                str = num2str(seq);
                fprintf(fid,'%s\n',str);
            end
            fclose(fid);


            file_name = ['PE_' num2str(num_PE) 'x' num2str(num_slices) '_etl' num2str(etl) '_f' num2str(accel_factor)  '_rep' num2str(rep) '_DCFill' num2str(central_region) '_pe2'];

            % write the second file
            fid = fopen(file_name, 'w');
            str = 't2 =';
            fprintf(fid,'%s\n',str);
            for i = 1:size(PE2,1)
                seq = PE2(i,:);
                str = num2str(seq);
                fprintf(fid,'%s\n',str);
            end
            fclose(fid);


        end
    end
end



