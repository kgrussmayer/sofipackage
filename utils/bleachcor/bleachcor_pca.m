function stack = bleachcor_pca(stack,settings)

    [DataCorrected,Background] = BgRemovalAutomatic(double(stack),2);
    stack = DataCorrected;

end
