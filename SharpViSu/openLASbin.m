function A = openLASbin (folder, filename)
        if all(filename(end-3:end) == '.bin')
            fileID = fopen([folder filename]);
            A = fread(fileID, 'uint32');
            A = reshape(A, 9, size(A,1)/9);
            A = A';
            fileID = fopen([folder filename]);
            B = fread(fileID, 'single');
            B = reshape(B, 9, size(B,1)/9);
            B = B';
            A(:,4:9) = B(:,4:9);
            fclose(fileID);
        end
end