1.运行代码：
打开demo_fatureselection
分别将其中的method改为“MRSR”，“RUFS”，“DPFS”，“AUFS”，“UDFS”，分别点击运行即可。(此处由于涉及版权问题，只提供MRSR的代码)。

2.代码说明：

（1）demo_featureselection.m
这个关于无监督选择的代码为一个自动化的框架，里面包含了多个方法的代码。
其中demo_featureselection为script文件，里面对于dataset,feature number,method,lamda等参数进行了设置，可以根据自己的需要进行修改。

（2）evalute.m
method为你要选择的无监督特征选择的方法，具体可以参见evalute.m函数。
evalute.m中使用switch case的方法给出了各种方法下的参数设置。

（3）evalute_num.m
evalute_num.m函数主要实现特征选择结果的评价，评价标准分别为classification，NMI，rebunduncy，clustering。


