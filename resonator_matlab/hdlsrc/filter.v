// -------------------------------------------------------------
//
// Module: filter
// Generated by MATLAB(R) 9.10 and Filter Design HDL Coder 3.1.9.
// Generated on: 2021-08-27 16:26:33
// -------------------------------------------------------------

// -------------------------------------------------------------
// HDL Code Generation Options:
//
// EDAScriptGeneration: off
// AddPipelineRegisters: on
// TargetLanguage: Verilog
// TestBenchStimulus: step ramp chirp 
// GenerateHDLTestBench: off

// Filter Specifications:
//
// 采样率  : 不适用(归一化频率)
// 响应   : Lowpass
// 设定   : Fp,Fst,Ap,Ast
// 阻带衰减 : 80 dB
// 阻带边缘 : 0.0008
// 通带边缘 : 2e-05
// 通带波纹 : 1 dB
// -------------------------------------------------------------

// -------------------------------------------------------------
// HDL Implementation    : Fully parallel
// Folding Factor        : 1
// -------------------------------------------------------------




`timescale 1 ns / 1 ns

module filter
               (
                clk,
                clk_enable,
                reset,
                filter_in,
                filter_out
                );

  input   clk; 
  input   clk_enable; 
  input   reset; 
  input   signed [15:0] filter_in; //sfix16_En15
  output  signed [31:0] filter_out; //sfix32_En27

////////////////////////////////////////////////////////////////
//Module Architecture: filter
////////////////////////////////////////////////////////////////
  // Local Functions
  // Type Definitions
  // Constants
  parameter signed [31:0] scaleconst1 = 32'h0001D390; //sfix32_En45
  parameter signed [31:0] coeff_b1_section1 = 32'h20000000; //sfix32_En29
  parameter signed [31:0] coeff_b2_section1 = 32'h40000000; //sfix32_En29
  parameter signed [31:0] coeff_b3_section1 = 32'h20000000; //sfix32_En29
  parameter signed [31:0] coeff_a2_section1 = 32'h8001E952; //sfix32_En30
  parameter signed [31:0] coeff_a3_section1 = 32'h3FFE16BD; //sfix32_En30
  parameter signed [31:0] scaleconst2 = 32'h7A50C1B7; //sfix32_En45
  parameter signed [31:0] coeff_b1_section2 = 32'h20000000; //sfix32_En29
  parameter signed [31:0] coeff_b2_section2 = 32'h20000000; //sfix32_En29
  parameter signed [31:0] coeff_b3_section2 = 32'h00000000; //sfix32_En29
  parameter signed [31:0] coeff_a2_section2 = 32'hC001E943; //sfix32_En30
  parameter signed [31:0] coeff_a3_section2 = 32'h00000000; //sfix32_En30
  // Signals
  reg  signed [15:0] input_register; // sfix16_En15
  wire signed [95:0] scale1; // sfix96_En72
  wire signed [47:0] mul_temp; // sfix48_En60
  wire signed [63:0] scaletypeconvert1; // sfix64_En40
  // Section 1 Signals 
  wire signed [97:0] a1sum1; // sfix98_En62
  wire signed [97:0] a2sum1; // sfix98_En62
  wire signed [97:0] b1sum1; // sfix98_En61
  wire signed [97:0] b2sum1; // sfix98_En61
  wire signed [63:0] typeconvert1; // sfix64_En32
  reg  signed [63:0] delay_section1 [0:1] ; // sfix64_En32
  wire signed [97:0] inputconv1; // sfix98_En62
  wire signed [95:0] a2mul1; // sfix96_En62
  wire signed [95:0] a3mul1; // sfix96_En62
  wire signed [95:0] b1mul1; // sfix96_En61
  wire signed [95:0] b2mul1; // sfix96_En61
  wire signed [95:0] b3mul1; // sfix96_En61
  wire signed [97:0] sub_cast; // sfix98_En62
  wire signed [97:0] sub_cast_1; // sfix98_En62
  wire signed [98:0] sub_temp; // sfix99_En62
  wire signed [97:0] sub_cast_2; // sfix98_En62
  wire signed [97:0] sub_cast_3; // sfix98_En62
  wire signed [98:0] sub_temp_1; // sfix99_En62
  wire signed [97:0] b1multypeconvert1; // sfix98_En61
  wire signed [97:0] add_cast; // sfix98_En61
  wire signed [97:0] add_cast_1; // sfix98_En61
  wire signed [98:0] add_temp; // sfix99_En61
  wire signed [97:0] add_cast_2; // sfix98_En61
  wire signed [97:0] add_cast_3; // sfix98_En61
  wire signed [98:0] add_temp_1; // sfix99_En61
  wire signed [63:0] section_result1; // sfix64_En27
  reg  signed [63:0] sos_pipeline1; // sfix64_En27
  wire signed [95:0] scale2; // sfix96_En72
  wire signed [63:0] scaletypeconvert2; // sfix64_En40
  //   -- Section 2 Signals 
  wire signed [97:0] a1sum2; // sfix98_En62
  wire signed [97:0] b1sum2; // sfix98_En61
  wire signed [63:0] a1sumtypeconvert2; // sfix64_En32
  reg  signed [63:0] delay_section2; // sfix64_En32
  wire signed [97:0] inputconv2; // sfix98_En62
  wire signed [95:0] a2mul2; // sfix96_En62
  wire signed [95:0] b1mul2; // sfix96_En61
  wire signed [95:0] b2mul2; // sfix96_En61
  wire signed [97:0] sub_cast_4; // sfix98_En62
  wire signed [97:0] sub_cast_5; // sfix98_En62
  wire signed [98:0] sub_temp_2; // sfix99_En62
  wire signed [97:0] b1multypeconvert2; // sfix98_En61
  wire signed [97:0] add_cast_4; // sfix98_En61
  wire signed [97:0] add_cast_5; // sfix98_En61
  wire signed [98:0] add_temp_2; // sfix99_En61
  wire signed [31:0] output_typeconvert; // sfix32_En27
  reg  signed [31:0] output_register; // sfix32_En27

  // Block Statements
  always @ (posedge clk or posedge reset)
    begin: input_reg_process
      if (reset == 1'b1) begin
        input_register <= 0;
      end
      else begin
        if (clk_enable == 1'b1) begin
          input_register <= filter_in;
        end
      end
    end // input_reg_process

  assign mul_temp = input_register * scaleconst1;
  assign scale1 = $signed({mul_temp[47:0], 12'b000000000000});

  assign scaletypeconvert1 = (scale1[95:0] + {scale1[32], {31{~scale1[32]}}})>>>32;

  //   ------------------ Section 1 ------------------

  assign typeconvert1 = (a1sum1[93:0] + {a1sum1[30], {29{~a1sum1[30]}}})>>>30;

  always @ (posedge clk or posedge reset)
    begin: delay_process_section1
      if (reset == 1'b1) begin
        delay_section1[0] <= 64'h0000000000000000;
        delay_section1[1] <= 64'h0000000000000000;
      end
      else begin
        if (clk_enable == 1'b1) begin
          delay_section1[1] <= delay_section1[0];
          delay_section1[0] <= typeconvert1;
        end
      end
    end // delay_process_section1

  assign inputconv1 = $signed({scaletypeconvert1[63:0], 22'b0000000000000000000000});

  assign a2mul1 = delay_section1[0] * coeff_a2_section1;

  assign a3mul1 = delay_section1[1] * coeff_a3_section1;

  assign b1mul1 = $signed({typeconvert1[63:0], 29'b00000000000000000000000000000});

  assign b2mul1 = $signed({delay_section1[0][63:0], 30'b000000000000000000000000000000});

  assign b3mul1 = $signed({delay_section1[1][63:0], 29'b00000000000000000000000000000});

  assign sub_cast = inputconv1;
  assign sub_cast_1 = $signed({{2{a2mul1[95]}}, a2mul1});
  assign sub_temp = sub_cast - sub_cast_1;
  assign a2sum1 = sub_temp[97:0];

  assign sub_cast_2 = a2sum1;
  assign sub_cast_3 = $signed({{2{a3mul1[95]}}, a3mul1});
  assign sub_temp_1 = sub_cast_2 - sub_cast_3;
  assign a1sum1 = sub_temp_1[97:0];

  assign b1multypeconvert1 = $signed({{2{b1mul1[95]}}, b1mul1});

  assign add_cast = b1multypeconvert1;
  assign add_cast_1 = $signed({{2{b2mul1[95]}}, b2mul1});
  assign add_temp = add_cast + add_cast_1;
  assign b2sum1 = add_temp[97:0];

  assign add_cast_2 = b2sum1;
  assign add_cast_3 = $signed({{2{b3mul1[95]}}, b3mul1});
  assign add_temp_1 = add_cast_2 + add_cast_3;
  assign b1sum1 = add_temp_1[97:0];

  assign section_result1 = (b1sum1[97:0] + {b1sum1[34], {33{~b1sum1[34]}}})>>>34;

  always @ (posedge clk or posedge reset)
    begin: sos_pipeline_process_section1
      if (reset == 1'b1) begin
        sos_pipeline1 <= 0;
      end
      else begin
        if (clk_enable == 1'b1) begin
          sos_pipeline1 <= section_result1;
        end
      end
    end // sos_pipeline_process_section1

  assign scale2 = sos_pipeline1 * scaleconst2;

  assign scaletypeconvert2 = (scale2[95:0] + {scale2[32], {31{~scale2[32]}}})>>>32;

  //   ------------------ Section 2 (First Order) ------------------

  assign a1sumtypeconvert2 = (a1sum2[93:0] + {a1sum2[30], {29{~a1sum2[30]}}})>>>30;

  always @ (posedge clk or posedge reset)
    begin: delay_process_section2
      if (reset == 1'b1) begin
        delay_section2 <= 0;
      end
      else begin
        if (clk_enable == 1'b1) begin
          delay_section2 <= a1sumtypeconvert2;
        end
      end
    end // delay_process_section2

  assign inputconv2 = $signed({scaletypeconvert2[63:0], 22'b0000000000000000000000});

  assign a2mul2 = delay_section2 * coeff_a2_section2;

  assign b1mul2 = $signed({a1sumtypeconvert2[63:0], 29'b00000000000000000000000000000});

  assign b2mul2 = $signed({delay_section2[63:0], 29'b00000000000000000000000000000});

  assign sub_cast_4 = inputconv2;
  assign sub_cast_5 = $signed({{2{a2mul2[95]}}, a2mul2});
  assign sub_temp_2 = sub_cast_4 - sub_cast_5;
  assign a1sum2 = sub_temp_2[97:0];

  assign b1multypeconvert2 = $signed({{2{b1mul2[95]}}, b1mul2});

  assign add_cast_4 = b1multypeconvert2;
  assign add_cast_5 = $signed({{2{b2mul2[95]}}, b2mul2});
  assign add_temp_2 = add_cast_4 + add_cast_5;
  assign b1sum2 = add_temp_2[97:0];

  assign output_typeconvert = (b1sum2[65:0] + {b1sum2[34], {33{~b1sum2[34]}}})>>>34;

  always @ (posedge clk or posedge reset)
    begin: Output_Register_process
      if (reset == 1'b1) begin
        output_register <= 0;
      end
      else begin
        if (clk_enable == 1'b1) begin
          output_register <= output_typeconvert;
        end
      end
    end // Output_Register_process

  // Assignment Statements
  assign filter_out = output_register;
endmodule  // filter
