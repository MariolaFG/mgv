#!/usr/bin/env python
# encoding: utf-8
import npyscreen
import curses
import sys
import re
import time

class ReadString(npyscreen.FixedText):
    def _print(self):
        string_to_print = self._get_string_to_print()
        if not string_to_print:
            return None
        string_to_print = string_to_print[self.begin_at:self.maximum_string_length+self.begin_at-self.left_margin]
        
        if sys.version_info[0] >= 3:
            string_to_print = self.display_value(self.value)[self.begin_at:self.maximum_string_length+self.begin_at-self.left_margin]
        else:
            # ensure unicode only here encoding here.
            dv = self.display_value(self.value)
            if isinstance(dv, bytes):
                dv = dv.decode(self.encoding, 'replace')
            string_to_print = dv[self.begin_at:self.maximum_string_length+self.begin_at-self.left_margin]
        
        column = 0
        place_in_string = 0
        if self.syntax_highlighting:
            self.update_highlighting(start=self.begin_at, end=self.maximum_string_length+self.begin_at-self.left_margin)
            while column <= (self.maximum_string_length - self.left_margin):
                if not string_to_print or place_in_string > len(string_to_print)-1:
                    break
                width_of_char_to_print = self.find_width_of_char(string_to_print[place_in_string])
                if column - 1 + width_of_char_to_print > self.maximum_string_length:
                    break 
                try:
                    highlight = self._highlightingdata[self.begin_at+place_in_string]
                except:
                    highlight = curses.A_NORMAL                
                self.parent.curses_pad.addstr(self.rely,self.relx+column+self.left_margin, 
                    self._print_unicode_char(string_to_print[place_in_string]), 
                    highlight
                    )
                column += self.find_width_of_char(string_to_print[place_in_string])
                place_in_string += 1
        else:
            if self.do_colors():
                if self.show_bold and self.color == 'DEFAULT':
                    color = self.parent.theme_manager.findPair(self, 'BOLD') | curses.A_BOLD
                elif self.show_bold:
                    color = self.parent.theme_manager.findPair(self, self.color) | curses.A_BOLD
                elif self.important:
                    color = self.parent.theme_manager.findPair(self, 'IMPORTANT') | curses.A_BOLD
                else:
                    color = self.parent.theme_manager.findPair(self)
            else:
                if self.important or self.show_bold:
                    color = curses.A_BOLD
                else:
                    color = curses.A_NORMAL

            pair = self.parent.theme_manager.get_pair_number("WHITE_BLACK")
            color= curses.color_pair(pair)
                   
            while column <= (self.maximum_string_length - self.left_margin):
                if not string_to_print or place_in_string > len(string_to_print)-1:
                    if self.highlight_whole_widget:
                        self.parent.curses_pad.addstr(self.rely, self.relx+column+self.left_margin, 
                                                      ' ',
                                                      color
                            )
                        column += width_of_char_to_print
                        place_in_string += 1
                        continue
                    else:
                        break
                        
                width_of_char_to_print = self.find_width_of_char(string_to_print[place_in_string])
                if column - 1 + width_of_char_to_print > self.maximum_string_length:
                    break

                if 1:
                    #str_chunks = re.split("([ACTG]| +)", self._print_unicode_char(string_to_print[place_in_string]))
                    str_chunks = re.split(SPLIT_TAG, string_to_print[place_in_string])
                    offset = 0
                    for chunk in str_chunks:
                        if chunk == '':
                            continue
                        elif chunk[0] == " ":
                            offset += len(chunk)
                            continue
                        elif chunk == "A":
                            pair = self.parent.theme_manager.get_pair_number("RED_BLACK")
                        elif chunk == "C":
                            pair = self.parent.theme_manager.get_pair_number("YELLOW_BLACK")
                        elif chunk == "T":
                            pair = self.parent.theme_manager.get_pair_number("CYAN_BLACK")
                        elif chunk == "G":
                            pair = self.parent.theme_manager.get_pair_number("GREEN_BLACK")
                        else:
                            pair = self.parent.theme_manager.get_pair_number("WHITE_BLACK")
                            
                        color= curses.color_pair(pair)

                        self.parent.curses_pad.addstr(self.rely, self.relx + column + self.left_margin + offset,
                                                      chunk,
                                                      color,
                                                      )
                        offset += len(chunk)
                else:
                    self.parent.curses_pad.addstr(self.rely, self.relx + column + self.left_margin,
                                                  self._print_unicode_char(string_to_print[place_in_string]),
                                                  color)

                column += width_of_char_to_print
                place_in_string += 1

#SPLIT_TAG = re.compile("([ACTG]| +)")
SPLIT_TAG = re.compile("(A+|C+|T+|G+| +|_+)")
class SeqView(npyscreen.MultiLine):
    _contained_widgets = ReadString

    def __init__(self, *args, **kargs):
        self.cursor_pos = 0
        npyscreen.MultiLine.__init__(self, *args, **kargs)
        self.max_length = max([len(v) for v in self.values])
        
    def update(self, clear=False):
        blank = curses.color_pair(self.parent.theme_manager.get_pair_number("WHITE_BLACK"))
        
        white = curses.color_pair(self.parent.theme_manager.get_pair_number("WHITE_BLACK"))
        red = curses.color_pair(self.parent.theme_manager.get_pair_number("RED_BLACK"))
        cyan = curses.color_pair(self.parent.theme_manager.get_pair_number("CYAN_BLACK"))
        green = curses.color_pair(self.parent.theme_manager.get_pair_number("GREEN_BLACK"))
        yellow = curses.color_pair(self.parent.theme_manager.get_pair_number("YELLOW_BLACK"))
        
        h_white = curses.color_pair(self.parent.theme_manager.get_pair_number("BLACK_WHITE"))
        h_red = curses.color_pair(self.parent.theme_manager.get_pair_number("BLACK_RED"))
        h_cyan = curses.color_pair(self.parent.theme_manager.get_pair_number("BLACK_CYAN"))
        h_green = curses.color_pair(self.parent.theme_manager.get_pair_number("BLACK_GREEN"))
        h_yellow = curses.color_pair(self.parent.theme_manager.get_pair_number("BLACK_YELLOW"))
        
        
        for rely, vl in enumerate(self.values[self.start_display_at:self.start_display_at+self.height-1]):
            trimmed_line = vl[self.cursor_pos : self.cursor_pos + self.width]
            str_chunks = re.split(SPLIT_TAG, trimmed_line)
            offset = 0
            hl = False
            if self.cursor_line == rely + self.start_display_at:
                hl = True
                self.parent.curses_pad.addstr(rely+1, 1, ">", yellow)
            else:
                self.parent.curses_pad.addstr(rely+1, 1, " ", blank)
                
            for chunk in str_chunks:
                if chunk == '':
                    continue
                elif chunk[0] == " ":
                    pair = blank
                elif chunk[0] == "_":
                    pair = yellow
                elif chunk == "A":
                    pair = red if not hl else h_red
                elif chunk == "C":
                    pair = yellow if not hl else h_yellow
                elif chunk == "T":
                    pair = cyan if not hl else h_cyan
                elif chunk == "G":
                    pair = green if not hl else h_green
                else:
                    pair = white if not hl else h_white

                self.parent.curses_pad.addstr(rely+1, 2+offset, chunk, pair)
                offset += len(chunk)
                
            if self.width > offset:
                self.parent.curses_pad.addstr(rely+1, 2+offset, " "*(self.width-offset), white)

    def display_value(self, vl):
        trimmed_line = vl[self.cursor_pos : self.cursor_pos + self.width]
        return npyscreen.MultiLine.display_value(self, trimmed_line)
        
    def move_right(self, x):
        self.cursor_pos = min(self.max_length, self.cursor_pos + 10)
        self.update(clear=True)

    def move_left(self, x):
        self.cursor_pos = max(0, self.cursor_pos -10)
        self.update(clear=True)
        
    def move_up(self, x):
        if self.cursor_line > 0:
            self.cursor_line -= 1
            
        if self.cursor_line < self.start_display_at:
            self.start_display_at = self.cursor_line
            
    def move_down(self, x):
        if self.cursor_line < len(self.values) - 1:
            self.cursor_line += 1
            
        if self.cursor_line > self.start_display_at + self.height:
            self.start_display_at += 1


        

        
class ReadViewApp(npyscreen.NPSApp):
    def __init__(self, values, title):
        self.values = values
        self.title = title
        npyscreen.NPSApp.__init__(self)
        
    def main(self):
        #npyscreen.setTheme(npyscreen.Themes.ColorfulTheme)
        F = npyscreen.Form(name = "Haplotype viewer: "+self.title)
        t3 = F.add(SeqView, always_show_cursor=True, values = self.values)

        def clean_exit(vl):
            t3.editing = False
            t3.how_exited = True
            F.ok_button.value = True
            
        t3.handlers.update({curses.KEY_RIGHT: t3.move_right,
                            curses.KEY_LEFT: t3.move_left,
                            curses.KEY_UP: t3.move_up,
                            curses.KEY_DOWN: t3.move_down,

                            
                            ord('q'): clean_exit})
        F.edit()

def view(values, title):
    App = ReadViewApp(values, title)
    App.run()
   
        
if __name__ == "__main__":
    import random
    S = ["A", "C", "T", "G"] 
    refseq = ''.join([S[random.randint(0, len(S)-1)] for i in xrange(5000)])

    temp_values = [refseq, refseq, refseq, refseq, refseq]

    R = "ACTG" + "-"*20
    
    for x in xrange(80):
        
        read = [R[random.randint(0, len(R)-1)] for i in xrange(random.randint(60, 120))]
        line = [" "] * len(refseq)
        start = random.randint(0, len(refseq)-len(read))
        line[start:start+len(read)] = read

        temp_values.append(''.join(line))
    view(temp_values, "")
    
