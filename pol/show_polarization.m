ikif true
    row=2048;  col=2448;
    fin=fopen('test/images/pol/90_bez_vzorku.raw','r');
    I=fread(fin,row*col,'uint8=>uint8'); 
    Z=reshape(I,col,row);
    Z=Z';
    % figure('Name',"Vstupni obrazek z RAWu")
    % k=imshow(Z);
    % title("Vstupni obrazek z RAWu")
end


J = demosaic(Z,"rggb");
%figure("Name","Demosaic")
%imshow(J)
%title("demosaic")

%% block pr
fun = @(block) mean(mean(block.data));

% fun = @(input) input;
proc = blockproc(Z, [2 2], fun);
figure("Name","Prumerovani a demosaic");
imshow(demosaic(uint8(proc),'rggb'),[]);
title("Prumerovani a demosaic")

%%
fun11 =  @(block) block.data(1,1);
fun12 =  @(block) block.data(1,2);
fun21 =  @(block) block.data(2,1);
fun22 =  @(block) block.data(2,2);

figure("Name","Polarizacni obrazky")
subplot(221)
proc1 = blockproc(Z, [2 2], fun11);
imshow(demosaic(uint8(proc1),'rggb'),[]);
title("90°")

subplot(222)
proc2 = blockproc(Z, [2 2], fun12);
imshow(demosaic(uint8(proc2),'rggb'),[]);
title("45°")

subplot(223)
proc3 = blockproc(Z, [2 2], fun21);
imshow(demosaic(uint8(proc3),'rggb'),[]);
title("135°")

subplot(224)
proc4 = blockproc(Z, [2 2], fun22);
imshow(demosaic(uint8(proc4),'rggb'),[]);
title("0°")
