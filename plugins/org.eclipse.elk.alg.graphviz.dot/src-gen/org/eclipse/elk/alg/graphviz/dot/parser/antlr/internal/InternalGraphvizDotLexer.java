/*******************************************************************************
 * Copyright (c) 2019 Kiel University and others.
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.

 * This Source Code may also be made available under the following Secondary
 * Licenses when the conditions for such availability set forth in the Eclipse
 * Public License v. 2.0 are satisfied: GPL-3.0 which is available at
 * https://www.gnu.org/licenses/gpl-3.0-standalone.html.
 * 
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.alg.graphviz.dot.parser.antlr.internal;

// Hack: Use our own Lexer superclass by means of import. 
// Currently there is no other way to specify the superclass for the lexer.
import org.eclipse.xtext.parser.antlr.Lexer;


import org.antlr.runtime.*;
import java.util.Stack;
import java.util.List;
import java.util.ArrayList;

@SuppressWarnings("all")
public class InternalGraphvizDotLexer extends Lexer {
    public static final int RULE_PREC_LINE=8;
    public static final int RULE_STRING=7;
    public static final int RULE_SL_COMMENT=10;
    public static final int T__19=19;
    public static final int T__15=15;
    public static final int T__16=16;
    public static final int T__17=17;
    public static final int T__18=18;
    public static final int T__13=13;
    public static final int T__14=14;
    public static final int EOF=-1;
    public static final int RULE_ID=4;
    public static final int RULE_WS=11;
    public static final int RULE_ANY_OTHER=12;
    public static final int T__26=26;
    public static final int T__27=27;
    public static final int T__28=28;
    public static final int RULE_INT=5;
    public static final int T__22=22;
    public static final int RULE_ML_COMMENT=9;
    public static final int T__23=23;
    public static final int RULE_FLOAT=6;
    public static final int T__24=24;
    public static final int T__25=25;
    public static final int T__20=20;
    public static final int T__21=21;

    // delegates
    // delegators

    public InternalGraphvizDotLexer() {;} 
    public InternalGraphvizDotLexer(CharStream input) {
        this(input, new RecognizerSharedState());
    }
    public InternalGraphvizDotLexer(CharStream input, RecognizerSharedState state) {
        super(input,state);

    }
    public String getGrammarFileName() { return "InternalGraphvizDot.g"; }

    // $ANTLR start "T__13"
    public final void mT__13() throws RecognitionException {
        try {
            int _type = T__13;
            int _channel = DEFAULT_TOKEN_CHANNEL;
            // InternalGraphvizDot.g:11:7: ( 'strict' )
            // InternalGraphvizDot.g:11:9: 'strict'
            {
            match("strict"); 


            }

            state.type = _type;
            state.channel = _channel;
        }
        finally {
        }
    }
    // $ANTLR end "T__13"

    // $ANTLR start "T__14"
    public final void mT__14() throws RecognitionException {
        try {
            int _type = T__14;
            int _channel = DEFAULT_TOKEN_CHANNEL;
            // InternalGraphvizDot.g:12:7: ( '{' )
            // InternalGraphvizDot.g:12:9: '{'
            {
            match('{'); 

            }

            state.type = _type;
            state.channel = _channel;
        }
        finally {
        }
    }
    // $ANTLR end "T__14"

    // $ANTLR start "T__15"
    public final void mT__15() throws RecognitionException {
        try {
            int _type = T__15;
            int _channel = DEFAULT_TOKEN_CHANNEL;
            // InternalGraphvizDot.g:13:7: ( '}' )
            // InternalGraphvizDot.g:13:9: '}'
            {
            match('}'); 

            }

            state.type = _type;
            state.channel = _channel;
        }
        finally {
        }
    }
    // $ANTLR end "T__15"

    // $ANTLR start "T__16"
    public final void mT__16() throws RecognitionException {
        try {
            int _type = T__16;
            int _channel = DEFAULT_TOKEN_CHANNEL;
            // InternalGraphvizDot.g:14:7: ( ';' )
            // InternalGraphvizDot.g:14:9: ';'
            {
            match(';'); 

            }

            state.type = _type;
            state.channel = _channel;
        }
        finally {
        }
    }
    // $ANTLR end "T__16"

    // $ANTLR start "T__17"
    public final void mT__17() throws RecognitionException {
        try {
            int _type = T__17;
            int _channel = DEFAULT_TOKEN_CHANNEL;
            // InternalGraphvizDot.g:15:7: ( '=' )
            // InternalGraphvizDot.g:15:9: '='
            {
            match('='); 

            }

            state.type = _type;
            state.channel = _channel;
        }
        finally {
        }
    }
    // $ANTLR end "T__17"

    // $ANTLR start "T__18"
    public final void mT__18() throws RecognitionException {
        try {
            int _type = T__18;
            int _channel = DEFAULT_TOKEN_CHANNEL;
            // InternalGraphvizDot.g:16:7: ( '[' )
            // InternalGraphvizDot.g:16:9: '['
            {
            match('['); 

            }

            state.type = _type;
            state.channel = _channel;
        }
        finally {
        }
    }
    // $ANTLR end "T__18"

    // $ANTLR start "T__19"
    public final void mT__19() throws RecognitionException {
        try {
            int _type = T__19;
            int _channel = DEFAULT_TOKEN_CHANNEL;
            // InternalGraphvizDot.g:17:7: ( ',' )
            // InternalGraphvizDot.g:17:9: ','
            {
            match(','); 

            }

            state.type = _type;
            state.channel = _channel;
        }
        finally {
        }
    }
    // $ANTLR end "T__19"

    // $ANTLR start "T__20"
    public final void mT__20() throws RecognitionException {
        try {
            int _type = T__20;
            int _channel = DEFAULT_TOKEN_CHANNEL;
            // InternalGraphvizDot.g:18:7: ( ']' )
            // InternalGraphvizDot.g:18:9: ']'
            {
            match(']'); 

            }

            state.type = _type;
            state.channel = _channel;
        }
        finally {
        }
    }
    // $ANTLR end "T__20"

    // $ANTLR start "T__21"
    public final void mT__21() throws RecognitionException {
        try {
            int _type = T__21;
            int _channel = DEFAULT_TOKEN_CHANNEL;
            // InternalGraphvizDot.g:19:7: ( ':' )
            // InternalGraphvizDot.g:19:9: ':'
            {
            match(':'); 

            }

            state.type = _type;
            state.channel = _channel;
        }
        finally {
        }
    }
    // $ANTLR end "T__21"

    // $ANTLR start "T__22"
    public final void mT__22() throws RecognitionException {
        try {
            int _type = T__22;
            int _channel = DEFAULT_TOKEN_CHANNEL;
            // InternalGraphvizDot.g:20:7: ( 'subgraph' )
            // InternalGraphvizDot.g:20:9: 'subgraph'
            {
            match("subgraph"); 


            }

            state.type = _type;
            state.channel = _channel;
        }
        finally {
        }
    }
    // $ANTLR end "T__22"

    // $ANTLR start "T__23"
    public final void mT__23() throws RecognitionException {
        try {
            int _type = T__23;
            int _channel = DEFAULT_TOKEN_CHANNEL;
            // InternalGraphvizDot.g:21:7: ( '->' )
            // InternalGraphvizDot.g:21:9: '->'
            {
            match("->"); 


            }

            state.type = _type;
            state.channel = _channel;
        }
        finally {
        }
    }
    // $ANTLR end "T__23"

    // $ANTLR start "T__24"
    public final void mT__24() throws RecognitionException {
        try {
            int _type = T__24;
            int _channel = DEFAULT_TOKEN_CHANNEL;
            // InternalGraphvizDot.g:22:7: ( '--' )
            // InternalGraphvizDot.g:22:9: '--'
            {
            match("--"); 


            }

            state.type = _type;
            state.channel = _channel;
        }
        finally {
        }
    }
    // $ANTLR end "T__24"

    // $ANTLR start "T__25"
    public final void mT__25() throws RecognitionException {
        try {
            int _type = T__25;
            int _channel = DEFAULT_TOKEN_CHANNEL;
            // InternalGraphvizDot.g:23:7: ( 'graph' )
            // InternalGraphvizDot.g:23:9: 'graph'
            {
            match("graph"); 


            }

            state.type = _type;
            state.channel = _channel;
        }
        finally {
        }
    }
    // $ANTLR end "T__25"

    // $ANTLR start "T__26"
    public final void mT__26() throws RecognitionException {
        try {
            int _type = T__26;
            int _channel = DEFAULT_TOKEN_CHANNEL;
            // InternalGraphvizDot.g:24:7: ( 'digraph' )
            // InternalGraphvizDot.g:24:9: 'digraph'
            {
            match("digraph"); 


            }

            state.type = _type;
            state.channel = _channel;
        }
        finally {
        }
    }
    // $ANTLR end "T__26"

    // $ANTLR start "T__27"
    public final void mT__27() throws RecognitionException {
        try {
            int _type = T__27;
            int _channel = DEFAULT_TOKEN_CHANNEL;
            // InternalGraphvizDot.g:25:7: ( 'node' )
            // InternalGraphvizDot.g:25:9: 'node'
            {
            match("node"); 


            }

            state.type = _type;
            state.channel = _channel;
        }
        finally {
        }
    }
    // $ANTLR end "T__27"

    // $ANTLR start "T__28"
    public final void mT__28() throws RecognitionException {
        try {
            int _type = T__28;
            int _channel = DEFAULT_TOKEN_CHANNEL;
            // InternalGraphvizDot.g:26:7: ( 'edge' )
            // InternalGraphvizDot.g:26:9: 'edge'
            {
            match("edge"); 


            }

            state.type = _type;
            state.channel = _channel;
        }
        finally {
        }
    }
    // $ANTLR end "T__28"

    // $ANTLR start "RULE_INT"
    public final void mRULE_INT() throws RecognitionException {
        try {
            int _type = RULE_INT;
            int _channel = DEFAULT_TOKEN_CHANNEL;
            // InternalGraphvizDot.g:1125:10: ( ( '-' )? ( '0' .. '9' )+ )
            // InternalGraphvizDot.g:1125:12: ( '-' )? ( '0' .. '9' )+
            {
            // InternalGraphvizDot.g:1125:12: ( '-' )?
            int alt1=2;
            int LA1_0 = input.LA(1);

            if ( (LA1_0=='-') ) {
                alt1=1;
            }
            switch (alt1) {
                case 1 :
                    // InternalGraphvizDot.g:1125:12: '-'
                    {
                    match('-'); 

                    }
                    break;

            }

            // InternalGraphvizDot.g:1125:17: ( '0' .. '9' )+
            int cnt2=0;
            loop2:
            do {
                int alt2=2;
                int LA2_0 = input.LA(1);

                if ( ((LA2_0>='0' && LA2_0<='9')) ) {
                    alt2=1;
                }


                switch (alt2) {
            	case 1 :
            	    // InternalGraphvizDot.g:1125:18: '0' .. '9'
            	    {
            	    matchRange('0','9'); 

            	    }
            	    break;

            	default :
            	    if ( cnt2 >= 1 ) break loop2;
                        EarlyExitException eee =
                            new EarlyExitException(2, input);
                        throw eee;
                }
                cnt2++;
            } while (true);


            }

            state.type = _type;
            state.channel = _channel;
        }
        finally {
        }
    }
    // $ANTLR end "RULE_INT"

    // $ANTLR start "RULE_ID"
    public final void mRULE_ID() throws RecognitionException {
        try {
            int _type = RULE_ID;
            int _channel = DEFAULT_TOKEN_CHANNEL;
            // InternalGraphvizDot.g:1127:9: ( ( 'a' .. 'z' | 'A' .. 'Z' | '_' ) ( 'a' .. 'z' | 'A' .. 'Z' | '_' | '0' .. '9' )* )
            // InternalGraphvizDot.g:1127:11: ( 'a' .. 'z' | 'A' .. 'Z' | '_' ) ( 'a' .. 'z' | 'A' .. 'Z' | '_' | '0' .. '9' )*
            {
            if ( (input.LA(1)>='A' && input.LA(1)<='Z')||input.LA(1)=='_'||(input.LA(1)>='a' && input.LA(1)<='z') ) {
                input.consume();

            }
            else {
                MismatchedSetException mse = new MismatchedSetException(null,input);
                recover(mse);
                throw mse;}

            // InternalGraphvizDot.g:1127:35: ( 'a' .. 'z' | 'A' .. 'Z' | '_' | '0' .. '9' )*
            loop3:
            do {
                int alt3=2;
                int LA3_0 = input.LA(1);

                if ( ((LA3_0>='0' && LA3_0<='9')||(LA3_0>='A' && LA3_0<='Z')||LA3_0=='_'||(LA3_0>='a' && LA3_0<='z')) ) {
                    alt3=1;
                }


                switch (alt3) {
            	case 1 :
            	    // InternalGraphvizDot.g:
            	    {
            	    if ( (input.LA(1)>='0' && input.LA(1)<='9')||(input.LA(1)>='A' && input.LA(1)<='Z')||input.LA(1)=='_'||(input.LA(1)>='a' && input.LA(1)<='z') ) {
            	        input.consume();

            	    }
            	    else {
            	        MismatchedSetException mse = new MismatchedSetException(null,input);
            	        recover(mse);
            	        throw mse;}


            	    }
            	    break;

            	default :
            	    break loop3;
                }
            } while (true);


            }

            state.type = _type;
            state.channel = _channel;
        }
        finally {
        }
    }
    // $ANTLR end "RULE_ID"

    // $ANTLR start "RULE_FLOAT"
    public final void mRULE_FLOAT() throws RecognitionException {
        try {
            int _type = RULE_FLOAT;
            int _channel = DEFAULT_TOKEN_CHANNEL;
            // InternalGraphvizDot.g:1129:12: ( ( '-' )? ( '0' .. '9' )* '.' ( '0' .. '9' )+ )
            // InternalGraphvizDot.g:1129:14: ( '-' )? ( '0' .. '9' )* '.' ( '0' .. '9' )+
            {
            // InternalGraphvizDot.g:1129:14: ( '-' )?
            int alt4=2;
            int LA4_0 = input.LA(1);

            if ( (LA4_0=='-') ) {
                alt4=1;
            }
            switch (alt4) {
                case 1 :
                    // InternalGraphvizDot.g:1129:14: '-'
                    {
                    match('-'); 

                    }
                    break;

            }

            // InternalGraphvizDot.g:1129:19: ( '0' .. '9' )*
            loop5:
            do {
                int alt5=2;
                int LA5_0 = input.LA(1);

                if ( ((LA5_0>='0' && LA5_0<='9')) ) {
                    alt5=1;
                }


                switch (alt5) {
            	case 1 :
            	    // InternalGraphvizDot.g:1129:20: '0' .. '9'
            	    {
            	    matchRange('0','9'); 

            	    }
            	    break;

            	default :
            	    break loop5;
                }
            } while (true);

            match('.'); 
            // InternalGraphvizDot.g:1129:35: ( '0' .. '9' )+
            int cnt6=0;
            loop6:
            do {
                int alt6=2;
                int LA6_0 = input.LA(1);

                if ( ((LA6_0>='0' && LA6_0<='9')) ) {
                    alt6=1;
                }


                switch (alt6) {
            	case 1 :
            	    // InternalGraphvizDot.g:1129:36: '0' .. '9'
            	    {
            	    matchRange('0','9'); 

            	    }
            	    break;

            	default :
            	    if ( cnt6 >= 1 ) break loop6;
                        EarlyExitException eee =
                            new EarlyExitException(6, input);
                        throw eee;
                }
                cnt6++;
            } while (true);


            }

            state.type = _type;
            state.channel = _channel;
        }
        finally {
        }
    }
    // $ANTLR end "RULE_FLOAT"

    // $ANTLR start "RULE_STRING"
    public final void mRULE_STRING() throws RecognitionException {
        try {
            int _type = RULE_STRING;
            int _channel = DEFAULT_TOKEN_CHANNEL;
            // InternalGraphvizDot.g:1131:13: ( '\"' ( '\\\\' '\"' | ~ ( '\"' ) )* '\"' )
            // InternalGraphvizDot.g:1131:15: '\"' ( '\\\\' '\"' | ~ ( '\"' ) )* '\"'
            {
            match('\"'); 
            // InternalGraphvizDot.g:1131:19: ( '\\\\' '\"' | ~ ( '\"' ) )*
            loop7:
            do {
                int alt7=3;
                int LA7_0 = input.LA(1);

                if ( (LA7_0=='\\') ) {
                    int LA7_2 = input.LA(2);

                    if ( (LA7_2=='\"') ) {
                        int LA7_4 = input.LA(3);

                        if ( ((LA7_4>='\u0000' && LA7_4<='\uFFFF')) ) {
                            alt7=1;
                        }

                        else {
                            alt7=2;
                        }

                    }
                    else if ( ((LA7_2>='\u0000' && LA7_2<='!')||(LA7_2>='#' && LA7_2<='\uFFFF')) ) {
                        alt7=2;
                    }


                }
                else if ( ((LA7_0>='\u0000' && LA7_0<='!')||(LA7_0>='#' && LA7_0<='[')||(LA7_0>=']' && LA7_0<='\uFFFF')) ) {
                    alt7=2;
                }


                switch (alt7) {
            	case 1 :
            	    // InternalGraphvizDot.g:1131:20: '\\\\' '\"'
            	    {
            	    match('\\'); 
            	    match('\"'); 

            	    }
            	    break;
            	case 2 :
            	    // InternalGraphvizDot.g:1131:29: ~ ( '\"' )
            	    {
            	    if ( (input.LA(1)>='\u0000' && input.LA(1)<='!')||(input.LA(1)>='#' && input.LA(1)<='\uFFFF') ) {
            	        input.consume();

            	    }
            	    else {
            	        MismatchedSetException mse = new MismatchedSetException(null,input);
            	        recover(mse);
            	        throw mse;}


            	    }
            	    break;

            	default :
            	    break loop7;
                }
            } while (true);

            match('\"'); 

            }

            state.type = _type;
            state.channel = _channel;
        }
        finally {
        }
    }
    // $ANTLR end "RULE_STRING"

    // $ANTLR start "RULE_PREC_LINE"
    public final void mRULE_PREC_LINE() throws RecognitionException {
        try {
            int _type = RULE_PREC_LINE;
            int _channel = DEFAULT_TOKEN_CHANNEL;
            // InternalGraphvizDot.g:1133:16: ( '#' (~ ( ( '\\n' | '\\r' ) ) )* ( ( '\\r' )? '\\n' )? )
            // InternalGraphvizDot.g:1133:18: '#' (~ ( ( '\\n' | '\\r' ) ) )* ( ( '\\r' )? '\\n' )?
            {
            match('#'); 
            // InternalGraphvizDot.g:1133:22: (~ ( ( '\\n' | '\\r' ) ) )*
            loop8:
            do {
                int alt8=2;
                int LA8_0 = input.LA(1);

                if ( ((LA8_0>='\u0000' && LA8_0<='\t')||(LA8_0>='\u000B' && LA8_0<='\f')||(LA8_0>='\u000E' && LA8_0<='\uFFFF')) ) {
                    alt8=1;
                }


                switch (alt8) {
            	case 1 :
            	    // InternalGraphvizDot.g:1133:22: ~ ( ( '\\n' | '\\r' ) )
            	    {
            	    if ( (input.LA(1)>='\u0000' && input.LA(1)<='\t')||(input.LA(1)>='\u000B' && input.LA(1)<='\f')||(input.LA(1)>='\u000E' && input.LA(1)<='\uFFFF') ) {
            	        input.consume();

            	    }
            	    else {
            	        MismatchedSetException mse = new MismatchedSetException(null,input);
            	        recover(mse);
            	        throw mse;}


            	    }
            	    break;

            	default :
            	    break loop8;
                }
            } while (true);

            // InternalGraphvizDot.g:1133:38: ( ( '\\r' )? '\\n' )?
            int alt10=2;
            int LA10_0 = input.LA(1);

            if ( (LA10_0=='\n'||LA10_0=='\r') ) {
                alt10=1;
            }
            switch (alt10) {
                case 1 :
                    // InternalGraphvizDot.g:1133:39: ( '\\r' )? '\\n'
                    {
                    // InternalGraphvizDot.g:1133:39: ( '\\r' )?
                    int alt9=2;
                    int LA9_0 = input.LA(1);

                    if ( (LA9_0=='\r') ) {
                        alt9=1;
                    }
                    switch (alt9) {
                        case 1 :
                            // InternalGraphvizDot.g:1133:39: '\\r'
                            {
                            match('\r'); 

                            }
                            break;

                    }

                    match('\n'); 

                    }
                    break;

            }


            }

            state.type = _type;
            state.channel = _channel;
        }
        finally {
        }
    }
    // $ANTLR end "RULE_PREC_LINE"

    // $ANTLR start "RULE_ML_COMMENT"
    public final void mRULE_ML_COMMENT() throws RecognitionException {
        try {
            int _type = RULE_ML_COMMENT;
            int _channel = DEFAULT_TOKEN_CHANNEL;
            // InternalGraphvizDot.g:1135:17: ( '/*' ( options {greedy=false; } : . )* '*/' )
            // InternalGraphvizDot.g:1135:19: '/*' ( options {greedy=false; } : . )* '*/'
            {
            match("/*"); 

            // InternalGraphvizDot.g:1135:24: ( options {greedy=false; } : . )*
            loop11:
            do {
                int alt11=2;
                int LA11_0 = input.LA(1);

                if ( (LA11_0=='*') ) {
                    int LA11_1 = input.LA(2);

                    if ( (LA11_1=='/') ) {
                        alt11=2;
                    }
                    else if ( ((LA11_1>='\u0000' && LA11_1<='.')||(LA11_1>='0' && LA11_1<='\uFFFF')) ) {
                        alt11=1;
                    }


                }
                else if ( ((LA11_0>='\u0000' && LA11_0<=')')||(LA11_0>='+' && LA11_0<='\uFFFF')) ) {
                    alt11=1;
                }


                switch (alt11) {
            	case 1 :
            	    // InternalGraphvizDot.g:1135:52: .
            	    {
            	    matchAny(); 

            	    }
            	    break;

            	default :
            	    break loop11;
                }
            } while (true);

            match("*/"); 


            }

            state.type = _type;
            state.channel = _channel;
        }
        finally {
        }
    }
    // $ANTLR end "RULE_ML_COMMENT"

    // $ANTLR start "RULE_SL_COMMENT"
    public final void mRULE_SL_COMMENT() throws RecognitionException {
        try {
            int _type = RULE_SL_COMMENT;
            int _channel = DEFAULT_TOKEN_CHANNEL;
            // InternalGraphvizDot.g:1137:17: ( '//' (~ ( ( '\\n' | '\\r' ) ) )* ( ( '\\r' )? '\\n' )? )
            // InternalGraphvizDot.g:1137:19: '//' (~ ( ( '\\n' | '\\r' ) ) )* ( ( '\\r' )? '\\n' )?
            {
            match("//"); 

            // InternalGraphvizDot.g:1137:24: (~ ( ( '\\n' | '\\r' ) ) )*
            loop12:
            do {
                int alt12=2;
                int LA12_0 = input.LA(1);

                if ( ((LA12_0>='\u0000' && LA12_0<='\t')||(LA12_0>='\u000B' && LA12_0<='\f')||(LA12_0>='\u000E' && LA12_0<='\uFFFF')) ) {
                    alt12=1;
                }


                switch (alt12) {
            	case 1 :
            	    // InternalGraphvizDot.g:1137:24: ~ ( ( '\\n' | '\\r' ) )
            	    {
            	    if ( (input.LA(1)>='\u0000' && input.LA(1)<='\t')||(input.LA(1)>='\u000B' && input.LA(1)<='\f')||(input.LA(1)>='\u000E' && input.LA(1)<='\uFFFF') ) {
            	        input.consume();

            	    }
            	    else {
            	        MismatchedSetException mse = new MismatchedSetException(null,input);
            	        recover(mse);
            	        throw mse;}


            	    }
            	    break;

            	default :
            	    break loop12;
                }
            } while (true);

            // InternalGraphvizDot.g:1137:40: ( ( '\\r' )? '\\n' )?
            int alt14=2;
            int LA14_0 = input.LA(1);

            if ( (LA14_0=='\n'||LA14_0=='\r') ) {
                alt14=1;
            }
            switch (alt14) {
                case 1 :
                    // InternalGraphvizDot.g:1137:41: ( '\\r' )? '\\n'
                    {
                    // InternalGraphvizDot.g:1137:41: ( '\\r' )?
                    int alt13=2;
                    int LA13_0 = input.LA(1);

                    if ( (LA13_0=='\r') ) {
                        alt13=1;
                    }
                    switch (alt13) {
                        case 1 :
                            // InternalGraphvizDot.g:1137:41: '\\r'
                            {
                            match('\r'); 

                            }
                            break;

                    }

                    match('\n'); 

                    }
                    break;

            }


            }

            state.type = _type;
            state.channel = _channel;
        }
        finally {
        }
    }
    // $ANTLR end "RULE_SL_COMMENT"

    // $ANTLR start "RULE_WS"
    public final void mRULE_WS() throws RecognitionException {
        try {
            int _type = RULE_WS;
            int _channel = DEFAULT_TOKEN_CHANNEL;
            // InternalGraphvizDot.g:1139:9: ( ( ' ' | '\\t' | '\\r' | '\\n' )+ )
            // InternalGraphvizDot.g:1139:11: ( ' ' | '\\t' | '\\r' | '\\n' )+
            {
            // InternalGraphvizDot.g:1139:11: ( ' ' | '\\t' | '\\r' | '\\n' )+
            int cnt15=0;
            loop15:
            do {
                int alt15=2;
                int LA15_0 = input.LA(1);

                if ( ((LA15_0>='\t' && LA15_0<='\n')||LA15_0=='\r'||LA15_0==' ') ) {
                    alt15=1;
                }


                switch (alt15) {
            	case 1 :
            	    // InternalGraphvizDot.g:
            	    {
            	    if ( (input.LA(1)>='\t' && input.LA(1)<='\n')||input.LA(1)=='\r'||input.LA(1)==' ' ) {
            	        input.consume();

            	    }
            	    else {
            	        MismatchedSetException mse = new MismatchedSetException(null,input);
            	        recover(mse);
            	        throw mse;}


            	    }
            	    break;

            	default :
            	    if ( cnt15 >= 1 ) break loop15;
                        EarlyExitException eee =
                            new EarlyExitException(15, input);
                        throw eee;
                }
                cnt15++;
            } while (true);


            }

            state.type = _type;
            state.channel = _channel;
        }
        finally {
        }
    }
    // $ANTLR end "RULE_WS"

    // $ANTLR start "RULE_ANY_OTHER"
    public final void mRULE_ANY_OTHER() throws RecognitionException {
        try {
            int _type = RULE_ANY_OTHER;
            int _channel = DEFAULT_TOKEN_CHANNEL;
            // InternalGraphvizDot.g:1141:16: ( . )
            // InternalGraphvizDot.g:1141:18: .
            {
            matchAny(); 

            }

            state.type = _type;
            state.channel = _channel;
        }
        finally {
        }
    }
    // $ANTLR end "RULE_ANY_OTHER"

    public void mTokens() throws RecognitionException {
        // InternalGraphvizDot.g:1:8: ( T__13 | T__14 | T__15 | T__16 | T__17 | T__18 | T__19 | T__20 | T__21 | T__22 | T__23 | T__24 | T__25 | T__26 | T__27 | T__28 | RULE_INT | RULE_ID | RULE_FLOAT | RULE_STRING | RULE_PREC_LINE | RULE_ML_COMMENT | RULE_SL_COMMENT | RULE_WS | RULE_ANY_OTHER )
        int alt16=25;
        alt16 = dfa16.predict(input);
        switch (alt16) {
            case 1 :
                // InternalGraphvizDot.g:1:10: T__13
                {
                mT__13(); 

                }
                break;
            case 2 :
                // InternalGraphvizDot.g:1:16: T__14
                {
                mT__14(); 

                }
                break;
            case 3 :
                // InternalGraphvizDot.g:1:22: T__15
                {
                mT__15(); 

                }
                break;
            case 4 :
                // InternalGraphvizDot.g:1:28: T__16
                {
                mT__16(); 

                }
                break;
            case 5 :
                // InternalGraphvizDot.g:1:34: T__17
                {
                mT__17(); 

                }
                break;
            case 6 :
                // InternalGraphvizDot.g:1:40: T__18
                {
                mT__18(); 

                }
                break;
            case 7 :
                // InternalGraphvizDot.g:1:46: T__19
                {
                mT__19(); 

                }
                break;
            case 8 :
                // InternalGraphvizDot.g:1:52: T__20
                {
                mT__20(); 

                }
                break;
            case 9 :
                // InternalGraphvizDot.g:1:58: T__21
                {
                mT__21(); 

                }
                break;
            case 10 :
                // InternalGraphvizDot.g:1:64: T__22
                {
                mT__22(); 

                }
                break;
            case 11 :
                // InternalGraphvizDot.g:1:70: T__23
                {
                mT__23(); 

                }
                break;
            case 12 :
                // InternalGraphvizDot.g:1:76: T__24
                {
                mT__24(); 

                }
                break;
            case 13 :
                // InternalGraphvizDot.g:1:82: T__25
                {
                mT__25(); 

                }
                break;
            case 14 :
                // InternalGraphvizDot.g:1:88: T__26
                {
                mT__26(); 

                }
                break;
            case 15 :
                // InternalGraphvizDot.g:1:94: T__27
                {
                mT__27(); 

                }
                break;
            case 16 :
                // InternalGraphvizDot.g:1:100: T__28
                {
                mT__28(); 

                }
                break;
            case 17 :
                // InternalGraphvizDot.g:1:106: RULE_INT
                {
                mRULE_INT(); 

                }
                break;
            case 18 :
                // InternalGraphvizDot.g:1:115: RULE_ID
                {
                mRULE_ID(); 

                }
                break;
            case 19 :
                // InternalGraphvizDot.g:1:123: RULE_FLOAT
                {
                mRULE_FLOAT(); 

                }
                break;
            case 20 :
                // InternalGraphvizDot.g:1:134: RULE_STRING
                {
                mRULE_STRING(); 

                }
                break;
            case 21 :
                // InternalGraphvizDot.g:1:146: RULE_PREC_LINE
                {
                mRULE_PREC_LINE(); 

                }
                break;
            case 22 :
                // InternalGraphvizDot.g:1:161: RULE_ML_COMMENT
                {
                mRULE_ML_COMMENT(); 

                }
                break;
            case 23 :
                // InternalGraphvizDot.g:1:177: RULE_SL_COMMENT
                {
                mRULE_SL_COMMENT(); 

                }
                break;
            case 24 :
                // InternalGraphvizDot.g:1:193: RULE_WS
                {
                mRULE_WS(); 

                }
                break;
            case 25 :
                // InternalGraphvizDot.g:1:201: RULE_ANY_OTHER
                {
                mRULE_ANY_OTHER(); 

                }
                break;

        }

    }


    protected DFA16 dfa16 = new DFA16(this);
    static final String DFA16_eotS =
        "\1\uffff\1\31\10\uffff\1\26\4\31\1\52\1\uffff\2\26\1\uffff\1\26\2\uffff\2\31\13\uffff\1\52\1\uffff\4\31\6\uffff\12\31\1\100\1\101\2\31\1\104\1\31\2\uffff\1\106\1\31\1\uffff\1\31\1\uffff\1\31\1\112\1\113\2\uffff";
    static final String DFA16_eofS =
        "\114\uffff";
    static final String DFA16_minS =
        "\1\0\1\164\10\uffff\1\55\1\162\1\151\1\157\1\144\1\56\1\uffff\1\60\1\0\1\uffff\1\52\2\uffff\1\162\1\142\13\uffff\1\56\1\uffff\1\141\1\147\1\144\1\147\6\uffff\1\151\1\147\1\160\1\162\2\145\1\143\1\162\1\150\1\141\2\60\1\164\1\141\1\60\1\160\2\uffff\1\60\1\160\1\uffff\1\150\1\uffff\1\150\2\60\2\uffff";
    static final String DFA16_maxS =
        "\1\uffff\1\165\10\uffff\1\76\1\162\1\151\1\157\1\144\1\71\1\uffff\1\71\1\uffff\1\uffff\1\57\2\uffff\1\162\1\142\13\uffff\1\71\1\uffff\1\141\1\147\1\144\1\147\6\uffff\1\151\1\147\1\160\1\162\2\145\1\143\1\162\1\150\1\141\2\172\1\164\1\141\1\172\1\160\2\uffff\1\172\1\160\1\uffff\1\150\1\uffff\1\150\2\172\2\uffff";
    static final String DFA16_acceptS =
        "\2\uffff\1\2\1\3\1\4\1\5\1\6\1\7\1\10\1\11\6\uffff\1\22\2\uffff\1\25\1\uffff\1\30\1\31\2\uffff\1\22\1\2\1\3\1\4\1\5\1\6\1\7\1\10\1\11\1\13\1\14\1\uffff\1\23\4\uffff\1\21\1\24\1\25\1\26\1\27\1\30\20\uffff\1\17\1\20\2\uffff\1\15\1\uffff\1\1\3\uffff\1\16\1\12";
    static final String DFA16_specialS =
        "\1\0\21\uffff\1\1\71\uffff}>";
    static final String[] DFA16_transitionS = {
            "\11\26\2\25\2\26\1\25\22\26\1\25\1\26\1\22\1\23\10\26\1\7\1\12\1\21\1\24\12\17\1\11\1\4\1\26\1\5\3\26\32\20\1\6\1\26\1\10\1\26\1\20\1\26\3\20\1\14\1\16\1\20\1\13\6\20\1\15\4\20\1\1\7\20\1\2\1\26\1\3\uff82\26",
            "\1\27\1\30",
            "",
            "",
            "",
            "",
            "",
            "",
            "",
            "",
            "\1\43\1\45\1\uffff\12\44\4\uffff\1\42",
            "\1\46",
            "\1\47",
            "\1\50",
            "\1\51",
            "\1\45\1\uffff\12\44",
            "",
            "\12\45",
            "\0\53",
            "",
            "\1\55\4\uffff\1\56",
            "",
            "",
            "\1\60",
            "\1\61",
            "",
            "",
            "",
            "",
            "",
            "",
            "",
            "",
            "",
            "",
            "",
            "\1\45\1\uffff\12\44",
            "",
            "\1\62",
            "\1\63",
            "\1\64",
            "\1\65",
            "",
            "",
            "",
            "",
            "",
            "",
            "\1\66",
            "\1\67",
            "\1\70",
            "\1\71",
            "\1\72",
            "\1\73",
            "\1\74",
            "\1\75",
            "\1\76",
            "\1\77",
            "\12\31\7\uffff\32\31\4\uffff\1\31\1\uffff\32\31",
            "\12\31\7\uffff\32\31\4\uffff\1\31\1\uffff\32\31",
            "\1\102",
            "\1\103",
            "\12\31\7\uffff\32\31\4\uffff\1\31\1\uffff\32\31",
            "\1\105",
            "",
            "",
            "\12\31\7\uffff\32\31\4\uffff\1\31\1\uffff\32\31",
            "\1\107",
            "",
            "\1\110",
            "",
            "\1\111",
            "\12\31\7\uffff\32\31\4\uffff\1\31\1\uffff\32\31",
            "\12\31\7\uffff\32\31\4\uffff\1\31\1\uffff\32\31",
            "",
            ""
    };

    static final short[] DFA16_eot = DFA.unpackEncodedString(DFA16_eotS);
    static final short[] DFA16_eof = DFA.unpackEncodedString(DFA16_eofS);
    static final char[] DFA16_min = DFA.unpackEncodedStringToUnsignedChars(DFA16_minS);
    static final char[] DFA16_max = DFA.unpackEncodedStringToUnsignedChars(DFA16_maxS);
    static final short[] DFA16_accept = DFA.unpackEncodedString(DFA16_acceptS);
    static final short[] DFA16_special = DFA.unpackEncodedString(DFA16_specialS);
    static final short[][] DFA16_transition;

    static {
        int numStates = DFA16_transitionS.length;
        DFA16_transition = new short[numStates][];
        for (int i=0; i<numStates; i++) {
            DFA16_transition[i] = DFA.unpackEncodedString(DFA16_transitionS[i]);
        }
    }

    class DFA16 extends DFA {

        public DFA16(BaseRecognizer recognizer) {
            this.recognizer = recognizer;
            this.decisionNumber = 16;
            this.eot = DFA16_eot;
            this.eof = DFA16_eof;
            this.min = DFA16_min;
            this.max = DFA16_max;
            this.accept = DFA16_accept;
            this.special = DFA16_special;
            this.transition = DFA16_transition;
        }
        public String getDescription() {
            return "1:1: Tokens : ( T__13 | T__14 | T__15 | T__16 | T__17 | T__18 | T__19 | T__20 | T__21 | T__22 | T__23 | T__24 | T__25 | T__26 | T__27 | T__28 | RULE_INT | RULE_ID | RULE_FLOAT | RULE_STRING | RULE_PREC_LINE | RULE_ML_COMMENT | RULE_SL_COMMENT | RULE_WS | RULE_ANY_OTHER );";
        }
        public int specialStateTransition(int s, IntStream _input) throws NoViableAltException {
            IntStream input = _input;
        	int _s = s;
            switch ( s ) {
                    case 0 : 
                        int LA16_0 = input.LA(1);

                        s = -1;
                        if ( (LA16_0=='s') ) {s = 1;}

                        else if ( (LA16_0=='{') ) {s = 2;}

                        else if ( (LA16_0=='}') ) {s = 3;}

                        else if ( (LA16_0==';') ) {s = 4;}

                        else if ( (LA16_0=='=') ) {s = 5;}

                        else if ( (LA16_0=='[') ) {s = 6;}

                        else if ( (LA16_0==',') ) {s = 7;}

                        else if ( (LA16_0==']') ) {s = 8;}

                        else if ( (LA16_0==':') ) {s = 9;}

                        else if ( (LA16_0=='-') ) {s = 10;}

                        else if ( (LA16_0=='g') ) {s = 11;}

                        else if ( (LA16_0=='d') ) {s = 12;}

                        else if ( (LA16_0=='n') ) {s = 13;}

                        else if ( (LA16_0=='e') ) {s = 14;}

                        else if ( ((LA16_0>='0' && LA16_0<='9')) ) {s = 15;}

                        else if ( ((LA16_0>='A' && LA16_0<='Z')||LA16_0=='_'||(LA16_0>='a' && LA16_0<='c')||LA16_0=='f'||(LA16_0>='h' && LA16_0<='m')||(LA16_0>='o' && LA16_0<='r')||(LA16_0>='t' && LA16_0<='z')) ) {s = 16;}

                        else if ( (LA16_0=='.') ) {s = 17;}

                        else if ( (LA16_0=='\"') ) {s = 18;}

                        else if ( (LA16_0=='#') ) {s = 19;}

                        else if ( (LA16_0=='/') ) {s = 20;}

                        else if ( ((LA16_0>='\t' && LA16_0<='\n')||LA16_0=='\r'||LA16_0==' ') ) {s = 21;}

                        else if ( ((LA16_0>='\u0000' && LA16_0<='\b')||(LA16_0>='\u000B' && LA16_0<='\f')||(LA16_0>='\u000E' && LA16_0<='\u001F')||LA16_0=='!'||(LA16_0>='$' && LA16_0<='+')||LA16_0=='<'||(LA16_0>='>' && LA16_0<='@')||LA16_0=='\\'||LA16_0=='^'||LA16_0=='`'||LA16_0=='|'||(LA16_0>='~' && LA16_0<='\uFFFF')) ) {s = 22;}

                        if ( s>=0 ) return s;
                        break;
                    case 1 : 
                        int LA16_18 = input.LA(1);

                        s = -1;
                        if ( ((LA16_18>='\u0000' && LA16_18<='\uFFFF')) ) {s = 43;}

                        else s = 22;

                        if ( s>=0 ) return s;
                        break;
            }
            NoViableAltException nvae =
                new NoViableAltException(getDescription(), 16, _s, input);
            error(nvae);
            throw nvae;
        }
    }
 

}