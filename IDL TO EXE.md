//保存为  .sav  文件

1 写好程序代码并编译运行无报错

2 在控制台运行命令 .FULL_RESET_SESSION

3 进行编译

4 在控制台运行命令  RESOLVE_ALL

5 在控制台运行命令 SAVE, /ROUTINES, FILENAME = 'myApp.sav'

//保存为 .exe 文件

6 接着运行 MAKE_RT, 'myApp', '/home/', savefile='myApp.sav', /OVERWRITE, /win32

7 进入编译文件夹打开 .ini文件 将第一行改为 Show=False




