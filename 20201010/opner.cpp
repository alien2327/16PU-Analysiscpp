#include <iostream>


int main() {
    //ﾌｧｲﾙを開くダイアログの作成 
    OpenFileDialog dlg = gcnew OpenFileDialog; 
    //ファイルフィルタ 
    dlg->Filter = "画像ﾌｧｲﾙ(*.bmp,*.jpg,*.png,*.tif,*.ico)|*.bmp;*.jpg;*.png;*.tif;*.ico"; 
    //ダイアログの表示 （Cancelボタンがクリックされた場合は何もしない）
    if (dlg->ShowDialog() == System::Windows::Forms::DialogResult::Cancel) return;
    //取得したファイル名 
    String FileName = dlg->FileName;
    return 0;
}
