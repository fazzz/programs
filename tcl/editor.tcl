#####################################################################################################
# メニューの作成
menu .m -type menubar
	. configure -menu .m
	.m add cascade -label "ファイル" -underline 0 -menu .m.m_fail 
	.m add cascade -label "計算" -underline 0 -menu .m.m_calc 
	.m add cascade -label "解析" -underline 0 -menu .m.m_anal 
	.m add cascade -label "表示" -underline 0 -menu .m.m_disp 
	.m add cascade -label "ヘルプ" -underline 0 -menu .m.m_help 

# ファイル・メニューの中身
menu .m.m_fail -tearoff no
.m.m_fail add command -label "Open" -under 0 -command "load_file"
.m.m_fail add separator
.m.m_fail add command -label "インプットファイルの作成" -under 0 -command "exec vis_mas_MD.exe"
.m.m_fail add separator
.m.m_fail add command -label "Exit" -under 0 -command exit
.m.m_fail add separator

set filename1 crd.in
set filename2 seq.in
set filename3 top.in
set filename4 clust.in
set mdinfilename md.in
set sequence "GLY GLY"

# 計算・メニューの中身
menu .m.m_calc -tearoff no
.m.m_calc add command -label "mas_MD" -command "exec mass_MD.exe -c $filename1 -r $filename4 -s $filename2 -t $filename3 -m $mdinfilename"
.m.m_calc add separator
.m.m_calc add command -label "vis_MD" -command "exec visdyn_pro7.exe -c coord.out -i isou.dat -m mol.dat"
.m.m_calc add separator
.m.m_calc add command -label "wat_MD" -command "exec md_test.exe -m md_w.in"
.m.m_calc add separator
.m.m_calc add command -label "pre_MD" -command "exec Vis_mass_MD.exe $sequence"
.m.m_calc add separator

# 解析・メニューの中身
menu .m.m_anal -tearoff no
.m.m_anal add command -label "make_energyprofile" -command "exec gawk -f pick_ene2.awk < thermo_dyn_properties.out"
.m.m_anal add separator
.m.m_anal add command -label "make_each_energyprofile" -command "exec gawk -f pick_ene.awk < thermo_dyn_properties.out"
.m.m_anal add separator
.m.m_anal add command -label "CalcDeltaHamilt" -command "exec ruby CalcDeltaHamilt.rb < thermo_dyn_properties.out > CalcDeltaHamilt.txt "
.m.m_anal add separator
.m.m_anal add command -label "rename" -command "exec mv *.txt *$mdinfilename.txt"
.m.m_anal add separator
.m.m_anal add command -label "excel" -command "exec excel"
.m.m_anal add separator

# 表示・メニューの中身
menu .m.m_disp -tearoff no
.m.m_disp add command -label "vis_MD" -command "exec visdyn_pro7.exe -c coord.out -i isou.dat -m mol.dat"
.m.m_disp add separator

# ヘルプ・メニューの中身
menu .m.m_help -tearoff no
###################################################################################################

###################################################################################################
# キャンパス
canvas .c0
###################################################################################################

###################################################################################################
#配列
frame .f_sequence

label .l_sequence -text "配列"
entry .e_sequence -textvariable sequence

#配置
#pack .l_sequence .e_sequence -in .f_sequence

###################################################################################################



###################################################################################################
#インプットファイル
frame .f_input
frame .f_input_ini
frame .f_input_seq
frame .f_input_top
frame .f_input_const

label .l_input_file -text "インプットファイル"

button .b_ini_coord -text "初期構造" -command "load_file_input1"
entry .e_ini_coord -textvariable filename1

set filename_seq ""
button .b_seq -text "sequence" -command "load_file_input2"
entry .e_seq -textvariable filename2

set filename_top ""
button .b_top -text "topology" -command "load_file_input3"
entry .e_top -textvariable filename3

set filename_clust ""
button .b_const -text "拘束条件" -command "load_file_input4"
entry .e_const -textvariable filename4

#配置
pack .b_ini_coord .e_ini_coord -in .f_input_ini
# -fil both
pack .b_seq .e_seq -in .f_input_seq 
# -fil both
pack .b_top .e_top -in .f_input_top 
# -fil both
pack .b_const .e_const -in .f_input_const 
# -fil both

pack .l_input_file .f_input_ini .f_input_seq .f_input_top .f_input_const -in .f_input
###################################################################################################

###################################################################################################
#計算条件
frame .f_input_calc
frame .f_input_mdin
frame .f_input_enns
frame .f_input_temp_initial
frame .f_input_temp_target
frame .f_input_inertia
frame .f_s_nvt
frame .f_dot_s_nvt
frame .f_input_numstep
frame .f_input_timestep
frame .f_input_outstep

set mdinfilename md.in

label .l_calc -text "計算条件"
label .l_mdin -text "ファイル名"
entry .e_mdin -textvariable mdinfilename

label .l_enns -text "統計集団"
radiobutton .rb_NVE -text "NVE" -variable enns -value 0 -command {
	destroy .l_temp_target .e_temp_target
	destroy .l_inertia .e_inertia
}
radiobutton .rb_NVT -text "NVT" -variable enns -value 1 -command {
	set temp_target  300
	label .l_temp_target -text "目標温度"
	entry .e_temp_target -textvariable temp_target
	if { $enns == 0} {
		entry .e_temp_target -state disabled
	}
	set inertia 1.0
	label .l_inertia -text "仮想系の緩和時間比"
	entry .e_inertia -textvariable inertia
	if { $enns == 0} {
		entry .e_inertia -state disabled
	}
	set s_nvt 1.0
	label .l_s_nvt -text "仮想系の変数"
	entry .e_s_nvt -textvariable s_NVT
	if { $enns == 0} {
		entry .e_s_nvt -state disabled
	}
	set dot_s_nvt 0.0
	label .l_dot_s_nvt -text "仮想系のの変数(速度)"
	entry .e_dot_s_nvt -textvariable dot_s_NVT
	if { $enns == 0} {
		entry .e_dot_s_nvt -state disabled
	}
	pack .l_temp_target .e_temp_target -in .f_input_temp_target
	pack .l_inertia .e_inertia -in .f_input_inertia
	pack .l_s_nvt .e_s_nvt -in .f_s_nvt
	pack .l_dot_s_nvt .e_dot_s_nvt -in .f_dot_s_nvt
}

set temp_initial 300
label .l_temp_initial -text "初期温度"
entry .e_temp_initial -textvariable temp_initial

#set temp_target  300
#label .l_temp_target -text "目標温度"
#entry .e_temp_target -textvariable temp_target
#if { $enns == 0} {
#	entry .e_temp_target -state disabled
#}
#
#
#set inertia 1.0
#label .l_inertia -text "仮想系の緩和時間比"
#entry .e_inertia -textvariable inertia
#if { $enns == 0} {
#	entry .e_inertia -state disabled
#}

set numstep 10000
label .l_numstep -text "ステップ数"
entry .e_numstep  -textvariable numstep

set timestep 0.001
label .l_timestep -text "時間幅"
entry .e_timestep -text "ps" -textvariable timestep

set outputstep 100
label .l_outstep -text "アウトプット幅"
entry .e_outstep -textvariable outputstep

button .b0 -text "確定" -command "mkin "

#配置
pack .l_mdin .e_mdin -in .f_input_mdin
pack .l_enns .rb_NVE .rb_NVT -in .f_input_enns
pack .l_temp_initial .e_temp_initial -in .f_input_temp_initial
#pack .l_temp_target .e_temp_target -in .f_input_temp_target
#pack .l_inertia .e_inertia -in .f_input_inertia
pack .l_numstep .e_numstep -in .f_input_numstep
pack .l_timestep .e_timestep -in .f_input_timestep
pack .l_outstep .e_outstep -in .f_input_outstep

pack .l_calc .f_input_mdin .f_input_enns .f_input_temp_initial .f_input_temp_target .f_input_inertia .f_s_nvt .f_dot_s_nvt .f_input_numstep .f_input_timestep .f_input_outstep  -in .f_input_calc
pack .b0 -side left -in .f_input_calc
###################################################################################################

###################################################################################################
#計算条件2
frame .f_input_calc2
frame .f_input_box
frame .f_input_cut

label .l_calc2 -text "計算条件2"

label .l_box -text "ボックスサイズ"
entry .e_box -textvariable box

label .l_cut -text "カットオフ半径"
entry .e_cut -textvariable cut

button .b1 -text "確定" -command "mkin2 "

#配置
pack .l_box .e_box -in .f_input_box
pack .l_cut .e_cut -in .f_input_cut

#pack .l_calc2 .f_input_box .f_input_cut -in .f_input_calc2
#pack .b1 -side left -in .f_input_calc2
###################################################################################################

###################################################################################################
# テキストウィジェット
frame .f_text
option add *font "{ＭＳ ゴシック} 10"
text .t_coord -xscrollcommand ".s_h set" -yscrollcommand ".s_v set" -wrap none
scrollbar .s_h -orient horizontal -command ".t_coord xview"
scrollbar .s_v -command ".t_coord yview"

set path_name ""
set num_flag 0
set line_count 0

label .l0  -anchor n -text "初期座標"

#配置
pack .l0 -in .f_text -anchor s
pack .s_v -side right -in .f_text
pack .s_h -side bottom -in .f_text
pack .t_coord -side bottom -anchor sw -in .f_text
###################################################################################################

###################################################################################################
#配置
#pack .f_text -side bottom
#pack .c0 -side left
pack .f_sequence -side left -fill both -padx 20
pack .f_input -side left -fill both -padx 40
pack .f_input_calc -side left -fill both -padx 20
pack .f_input_calc2 -side left -fill both -padx 20
###################################################################################################

###################################################################################################
proc load_file {} {
    global path_name num_flag line_count
    set filename [tk_getOpenFile -initialdir $path_name \
                                 -filetypes {{{InputFiles} {.pdb}}}]
    if [string compare $filename ""] {
        set path_name [file dirname $filename]
        # テキストをクリア
        .t_coord delete 1.0 end
        # ファイルの読み込み
        set f [open $filename]
        set line_count 1
        while {[gets $f line] >= 0} {
            .t_coord insert end "$line\n"
            incr line_count
        }
        close $f
        if $num_flag {
            insert_number
        }
    }
}
###################################################################################################

###################################################################################################
proc load_file_input1 {} {
    global filename1
    set filename1 [tk_getOpenFile -filetypes {{{InputFiles} {.in}}}]
}
###################################################################################################

###################################################################################################
proc load_file_input2 {} {
    global filename2
    set filename2 [tk_getOpenFile -filetypes {{{InputFiles} {.seq .in}}}]
}
###################################################################################################

###################################################################################################
proc load_file_input3 {} {
    global filename3
    set filename3 [tk_getOpenFile -filetypes {{{InputFiles} {.top .in}}}]
}
###################################################################################################

###################################################################################################
proc load_file_input4 {} {
    global filename4
    set filename4 [tk_getOpenFile -filetypes {{{InputFiles} {.clut .in}}}]
}
###################################################################################################

###################################################################################################
proc mkin {} {
	global enns temp_initial temp_target numstep timestep outputstep mdinfilename inertia s_nvt dot_s_nvt filename1 filename4 filename2 filename3
	set f [open $mdinfilename w+]
	puts $f $enns
	puts $f 1
	puts $f $temp_target
	puts $f $temp_initial
	puts $f $timestep
	puts $f $numstep
	puts $f 1
	puts $f 20
	puts $f $outputstep
	puts $f 1
	puts $f 0
	puts $f 0
	puts $f 0.0
	puts $f $inertia
	puts $f $s_nvt
	close $f

	set mdinfilename2 $mdinfilename.txt

	set f2 [open $mdinfilename2 w+]
	puts $f2 ---------------
	# 統計集団
	puts $f2 統計集団
	if {$enns == 0} {
		puts $f2 NVE
	} else {
		puts $f2 NVT
	}
	# タイムステップ、ステップ数、計算時間
	puts $f2 タイムステップ
	puts $f2 $timestep
	puts $f2 ステップ数
	puts $f2 $numstep
	# 初期温度
	puts $f2 初期温度
	puts $f2 $temp_initial
	if {$enns == 1} {
		# 目標温度
		puts $f2 目標温度
		puts $f2 $temp_target
		# 仮想慣性
		puts $f2 仮想慣性
		puts $f2 $inertia
		puts $f2 仮想変数
		puts $f2 $s_nvt
		puts $f2 仮想変数の速度
		puts $f2 $dot_s_nvt
	}
	puts $f2 ---------------
	puts $f2 初期構造ファイル
	puts $f2 $filename1
	puts $f2 剛体情報ファイル
	puts $f2 $filename4
	puts $f2 力場情報ファイル
	puts $f2 $filename2
	puts $f2 配列情報ファイル
	puts $f2 $filename3
	puts $f2 計算条件ファイル
	puts $f2 $mdinfilename
	puts $f2 ---------------
	close $f2
}
###################################################################################################

###################################################################################################
proc mkin2 {} {
	global box cut 
	set f [open "mdin2" w+]
	puts $f $box
	puts $f $cut
	close $f
}
###################################################################################################

