/*******************************************************************************
 * Copyright (c) 2019 Kiel University and others.
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 * 
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.graph.text.ide.contentassist.antlr.internal;

// Hack: Use our own Lexer superclass by means of import. 
// Currently there is no other way to specify the superclass for the lexer.
import org.eclipse.xtext.ide.editor.contentassist.antlr.internal.Lexer;


import org.antlr.runtime.*;
import java.util.Stack;
import java.util.List;
import java.util.ArrayList;

@SuppressWarnings("all")
public class InternalElkGraphLexer extends Lexer {
    public static final int RULE_STRING=4;
    public static final int RULE_SL_COMMENT=10;
    public static final int T__19=19;
    public static final int T__15=15;
    public static final int T__37=37;
    public static final int T__16=16;
    public static final int T__38=38;
    public static final int T__17=17;
    public static final int T__18=18;
    public static final int T__33=33;
    public static final int T__34=34;
    public static final int T__13=13;
    public static final int T__35=35;
    public static final int T__14=14;
    public static final int T__36=36;
    public static final int EOF=-1;
    public static final int T__30=30;
    public static final int T__31=31;
    public static final int T__32=32;
    public static final int RULE_ID=7;
    public static final int RULE_WS=11;
    public static final int RULE_SIGNED_INT=5;
    public static final int RULE_ANY_OTHER=12;
    public static final int T__26=26;
    public static final int T__27=27;
    public static final int T__28=28;
    public static final int RULE_INT=8;
    public static final int T__29=29;
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

    public InternalElkGraphLexer() {;} 
    public InternalElkGraphLexer(CharStream input) {
        this(input, new RecognizerSharedState());
    }
    public InternalElkGraphLexer(CharStream input, RecognizerSharedState state) {
        super(input,state);

    }
    public String getGrammarFileName() { return "InternalElkGraph.g"; }

    // $ANTLR start "T__13"
    public final void mT__13() throws RecognitionException {
        try {
            int _type = T__13;
            int _channel = DEFAULT_TOKEN_CHANNEL;
            // InternalElkGraph.g:11:7: ( 'null' )
            // InternalElkGraph.g:11:9: 'null'
            {
            match("null"); 


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
            // InternalElkGraph.g:12:7: ( 'true' )
            // InternalElkGraph.g:12:9: 'true'
            {
            match("true"); 


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
            // InternalElkGraph.g:13:7: ( 'false' )
            // InternalElkGraph.g:13:9: 'false'
            {
            match("false"); 


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
            // InternalElkGraph.g:14:7: ( 'graph' )
            // InternalElkGraph.g:14:9: 'graph'
            {
            match("graph"); 


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
            // InternalElkGraph.g:15:7: ( 'node' )
            // InternalElkGraph.g:15:9: 'node'
            {
            match("node"); 


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
            // InternalElkGraph.g:16:7: ( '{' )
            // InternalElkGraph.g:16:9: '{'
            {
            match('{'); 

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
            // InternalElkGraph.g:17:7: ( '}' )
            // InternalElkGraph.g:17:9: '}'
            {
            match('}'); 

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
            // InternalElkGraph.g:18:7: ( 'label' )
            // InternalElkGraph.g:18:9: 'label'
            {
            match("label"); 


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
            // InternalElkGraph.g:19:7: ( ':' )
            // InternalElkGraph.g:19:9: ':'
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
            // InternalElkGraph.g:20:7: ( 'port' )
            // InternalElkGraph.g:20:9: 'port'
            {
            match("port"); 


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
            // InternalElkGraph.g:21:7: ( 'layout' )
            // InternalElkGraph.g:21:9: 'layout'
            {
            match("layout"); 


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
            // InternalElkGraph.g:22:7: ( '[' )
            // InternalElkGraph.g:22:9: '['
            {
            match('['); 

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
            // InternalElkGraph.g:23:7: ( ']' )
            // InternalElkGraph.g:23:9: ']'
            {
            match(']'); 

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
            // InternalElkGraph.g:24:7: ( 'position' )
            // InternalElkGraph.g:24:9: 'position'
            {
            match("position"); 


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
            // InternalElkGraph.g:25:7: ( ',' )
            // InternalElkGraph.g:25:9: ','
            {
            match(','); 

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
            // InternalElkGraph.g:26:7: ( 'size' )
            // InternalElkGraph.g:26:9: 'size'
            {
            match("size"); 


            }

            state.type = _type;
            state.channel = _channel;
        }
        finally {
        }
    }
    // $ANTLR end "T__28"

    // $ANTLR start "T__29"
    public final void mT__29() throws RecognitionException {
        try {
            int _type = T__29;
            int _channel = DEFAULT_TOKEN_CHANNEL;
            // InternalElkGraph.g:27:7: ( 'edge' )
            // InternalElkGraph.g:27:9: 'edge'
            {
            match("edge"); 


            }

            state.type = _type;
            state.channel = _channel;
        }
        finally {
        }
    }
    // $ANTLR end "T__29"

    // $ANTLR start "T__30"
    public final void mT__30() throws RecognitionException {
        try {
            int _type = T__30;
            int _channel = DEFAULT_TOKEN_CHANNEL;
            // InternalElkGraph.g:28:7: ( '->' )
            // InternalElkGraph.g:28:9: '->'
            {
            match("->"); 


            }

            state.type = _type;
            state.channel = _channel;
        }
        finally {
        }
    }
    // $ANTLR end "T__30"

    // $ANTLR start "T__31"
    public final void mT__31() throws RecognitionException {
        try {
            int _type = T__31;
            int _channel = DEFAULT_TOKEN_CHANNEL;
            // InternalElkGraph.g:29:7: ( 'incoming' )
            // InternalElkGraph.g:29:9: 'incoming'
            {
            match("incoming"); 


            }

            state.type = _type;
            state.channel = _channel;
        }
        finally {
        }
    }
    // $ANTLR end "T__31"

    // $ANTLR start "T__32"
    public final void mT__32() throws RecognitionException {
        try {
            int _type = T__32;
            int _channel = DEFAULT_TOKEN_CHANNEL;
            // InternalElkGraph.g:30:7: ( 'outgoing' )
            // InternalElkGraph.g:30:9: 'outgoing'
            {
            match("outgoing"); 


            }

            state.type = _type;
            state.channel = _channel;
        }
        finally {
        }
    }
    // $ANTLR end "T__32"

    // $ANTLR start "T__33"
    public final void mT__33() throws RecognitionException {
        try {
            int _type = T__33;
            int _channel = DEFAULT_TOKEN_CHANNEL;
            // InternalElkGraph.g:31:7: ( 'start' )
            // InternalElkGraph.g:31:9: 'start'
            {
            match("start"); 


            }

            state.type = _type;
            state.channel = _channel;
        }
        finally {
        }
    }
    // $ANTLR end "T__33"

    // $ANTLR start "T__34"
    public final void mT__34() throws RecognitionException {
        try {
            int _type = T__34;
            int _channel = DEFAULT_TOKEN_CHANNEL;
            // InternalElkGraph.g:32:7: ( 'end' )
            // InternalElkGraph.g:32:9: 'end'
            {
            match("end"); 


            }

            state.type = _type;
            state.channel = _channel;
        }
        finally {
        }
    }
    // $ANTLR end "T__34"

    // $ANTLR start "T__35"
    public final void mT__35() throws RecognitionException {
        try {
            int _type = T__35;
            int _channel = DEFAULT_TOKEN_CHANNEL;
            // InternalElkGraph.g:33:7: ( 'bends' )
            // InternalElkGraph.g:33:9: 'bends'
            {
            match("bends"); 


            }

            state.type = _type;
            state.channel = _channel;
        }
        finally {
        }
    }
    // $ANTLR end "T__35"

    // $ANTLR start "T__36"
    public final void mT__36() throws RecognitionException {
        try {
            int _type = T__36;
            int _channel = DEFAULT_TOKEN_CHANNEL;
            // InternalElkGraph.g:34:7: ( '|' )
            // InternalElkGraph.g:34:9: '|'
            {
            match('|'); 

            }

            state.type = _type;
            state.channel = _channel;
        }
        finally {
        }
    }
    // $ANTLR end "T__36"

    // $ANTLR start "T__37"
    public final void mT__37() throws RecognitionException {
        try {
            int _type = T__37;
            int _channel = DEFAULT_TOKEN_CHANNEL;
            // InternalElkGraph.g:35:7: ( 'section' )
            // InternalElkGraph.g:35:9: 'section'
            {
            match("section"); 


            }

            state.type = _type;
            state.channel = _channel;
        }
        finally {
        }
    }
    // $ANTLR end "T__37"

    // $ANTLR start "T__38"
    public final void mT__38() throws RecognitionException {
        try {
            int _type = T__38;
            int _channel = DEFAULT_TOKEN_CHANNEL;
            // InternalElkGraph.g:36:7: ( '.' )
            // InternalElkGraph.g:36:9: '.'
            {
            match('.'); 

            }

            state.type = _type;
            state.channel = _channel;
        }
        finally {
        }
    }
    // $ANTLR end "T__38"

    // $ANTLR start "RULE_SIGNED_INT"
    public final void mRULE_SIGNED_INT() throws RecognitionException {
        try {
            int _type = RULE_SIGNED_INT;
            int _channel = DEFAULT_TOKEN_CHANNEL;
            // InternalElkGraph.g:5987:17: ( ( '+' | '-' )? RULE_INT )
            // InternalElkGraph.g:5987:19: ( '+' | '-' )? RULE_INT
            {
            // InternalElkGraph.g:5987:19: ( '+' | '-' )?
            int alt1=2;
            int LA1_0 = input.LA(1);

            if ( (LA1_0=='+'||LA1_0=='-') ) {
                alt1=1;
            }
            switch (alt1) {
                case 1 :
                    // InternalElkGraph.g:
                    {
                    if ( input.LA(1)=='+'||input.LA(1)=='-' ) {
                        input.consume();

                    }
                    else {
                        MismatchedSetException mse = new MismatchedSetException(null,input);
                        recover(mse);
                        throw mse;}


                    }
                    break;

            }

            mRULE_INT(); 

            }

            state.type = _type;
            state.channel = _channel;
        }
        finally {
        }
    }
    // $ANTLR end "RULE_SIGNED_INT"

    // $ANTLR start "RULE_FLOAT"
    public final void mRULE_FLOAT() throws RecognitionException {
        try {
            int _type = RULE_FLOAT;
            int _channel = DEFAULT_TOKEN_CHANNEL;
            // InternalElkGraph.g:5989:12: ( ( '+' | '-' )? ( RULE_INT '.' RULE_INT | RULE_INT ( '.' RULE_INT )? ( 'e' | 'E' ) ( '+' | '-' )? RULE_INT ) )
            // InternalElkGraph.g:5989:14: ( '+' | '-' )? ( RULE_INT '.' RULE_INT | RULE_INT ( '.' RULE_INT )? ( 'e' | 'E' ) ( '+' | '-' )? RULE_INT )
            {
            // InternalElkGraph.g:5989:14: ( '+' | '-' )?
            int alt2=2;
            int LA2_0 = input.LA(1);

            if ( (LA2_0=='+'||LA2_0=='-') ) {
                alt2=1;
            }
            switch (alt2) {
                case 1 :
                    // InternalElkGraph.g:
                    {
                    if ( input.LA(1)=='+'||input.LA(1)=='-' ) {
                        input.consume();

                    }
                    else {
                        MismatchedSetException mse = new MismatchedSetException(null,input);
                        recover(mse);
                        throw mse;}


                    }
                    break;

            }

            // InternalElkGraph.g:5989:25: ( RULE_INT '.' RULE_INT | RULE_INT ( '.' RULE_INT )? ( 'e' | 'E' ) ( '+' | '-' )? RULE_INT )
            int alt5=2;
            alt5 = dfa5.predict(input);
            switch (alt5) {
                case 1 :
                    // InternalElkGraph.g:5989:26: RULE_INT '.' RULE_INT
                    {
                    mRULE_INT(); 
                    match('.'); 
                    mRULE_INT(); 

                    }
                    break;
                case 2 :
                    // InternalElkGraph.g:5989:48: RULE_INT ( '.' RULE_INT )? ( 'e' | 'E' ) ( '+' | '-' )? RULE_INT
                    {
                    mRULE_INT(); 
                    // InternalElkGraph.g:5989:57: ( '.' RULE_INT )?
                    int alt3=2;
                    int LA3_0 = input.LA(1);

                    if ( (LA3_0=='.') ) {
                        alt3=1;
                    }
                    switch (alt3) {
                        case 1 :
                            // InternalElkGraph.g:5989:58: '.' RULE_INT
                            {
                            match('.'); 
                            mRULE_INT(); 

                            }
                            break;

                    }

                    if ( input.LA(1)=='E'||input.LA(1)=='e' ) {
                        input.consume();

                    }
                    else {
                        MismatchedSetException mse = new MismatchedSetException(null,input);
                        recover(mse);
                        throw mse;}

                    // InternalElkGraph.g:5989:83: ( '+' | '-' )?
                    int alt4=2;
                    int LA4_0 = input.LA(1);

                    if ( (LA4_0=='+'||LA4_0=='-') ) {
                        alt4=1;
                    }
                    switch (alt4) {
                        case 1 :
                            // InternalElkGraph.g:
                            {
                            if ( input.LA(1)=='+'||input.LA(1)=='-' ) {
                                input.consume();

                            }
                            else {
                                MismatchedSetException mse = new MismatchedSetException(null,input);
                                recover(mse);
                                throw mse;}


                            }
                            break;

                    }

                    mRULE_INT(); 

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
    // $ANTLR end "RULE_FLOAT"

    // $ANTLR start "RULE_ID"
    public final void mRULE_ID() throws RecognitionException {
        try {
            int _type = RULE_ID;
            int _channel = DEFAULT_TOKEN_CHANNEL;
            // InternalElkGraph.g:5991:9: ( ( '^' )? ( 'a' .. 'z' | 'A' .. 'Z' | '_' ) ( 'a' .. 'z' | 'A' .. 'Z' | '_' | '0' .. '9' )* )
            // InternalElkGraph.g:5991:11: ( '^' )? ( 'a' .. 'z' | 'A' .. 'Z' | '_' ) ( 'a' .. 'z' | 'A' .. 'Z' | '_' | '0' .. '9' )*
            {
            // InternalElkGraph.g:5991:11: ( '^' )?
            int alt6=2;
            int LA6_0 = input.LA(1);

            if ( (LA6_0=='^') ) {
                alt6=1;
            }
            switch (alt6) {
                case 1 :
                    // InternalElkGraph.g:5991:11: '^'
                    {
                    match('^'); 

                    }
                    break;

            }

            if ( (input.LA(1)>='A' && input.LA(1)<='Z')||input.LA(1)=='_'||(input.LA(1)>='a' && input.LA(1)<='z') ) {
                input.consume();

            }
            else {
                MismatchedSetException mse = new MismatchedSetException(null,input);
                recover(mse);
                throw mse;}

            // InternalElkGraph.g:5991:40: ( 'a' .. 'z' | 'A' .. 'Z' | '_' | '0' .. '9' )*
            loop7:
            do {
                int alt7=2;
                int LA7_0 = input.LA(1);

                if ( ((LA7_0>='0' && LA7_0<='9')||(LA7_0>='A' && LA7_0<='Z')||LA7_0=='_'||(LA7_0>='a' && LA7_0<='z')) ) {
                    alt7=1;
                }


                switch (alt7) {
            	case 1 :
            	    // InternalElkGraph.g:
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
            	    break loop7;
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

    // $ANTLR start "RULE_INT"
    public final void mRULE_INT() throws RecognitionException {
        try {
            // InternalElkGraph.g:5993:19: ( ( '0' .. '9' )+ )
            // InternalElkGraph.g:5993:21: ( '0' .. '9' )+
            {
            // InternalElkGraph.g:5993:21: ( '0' .. '9' )+
            int cnt8=0;
            loop8:
            do {
                int alt8=2;
                int LA8_0 = input.LA(1);

                if ( ((LA8_0>='0' && LA8_0<='9')) ) {
                    alt8=1;
                }


                switch (alt8) {
            	case 1 :
            	    // InternalElkGraph.g:5993:22: '0' .. '9'
            	    {
            	    matchRange('0','9'); 

            	    }
            	    break;

            	default :
            	    if ( cnt8 >= 1 ) break loop8;
                        EarlyExitException eee =
                            new EarlyExitException(8, input);
                        throw eee;
                }
                cnt8++;
            } while (true);


            }

        }
        finally {
        }
    }
    // $ANTLR end "RULE_INT"

    // $ANTLR start "RULE_STRING"
    public final void mRULE_STRING() throws RecognitionException {
        try {
            int _type = RULE_STRING;
            int _channel = DEFAULT_TOKEN_CHANNEL;
            // InternalElkGraph.g:5995:13: ( ( '\"' ( '\\\\' . | ~ ( ( '\\\\' | '\"' ) ) )* '\"' | '\\'' ( '\\\\' . | ~ ( ( '\\\\' | '\\'' ) ) )* '\\'' ) )
            // InternalElkGraph.g:5995:15: ( '\"' ( '\\\\' . | ~ ( ( '\\\\' | '\"' ) ) )* '\"' | '\\'' ( '\\\\' . | ~ ( ( '\\\\' | '\\'' ) ) )* '\\'' )
            {
            // InternalElkGraph.g:5995:15: ( '\"' ( '\\\\' . | ~ ( ( '\\\\' | '\"' ) ) )* '\"' | '\\'' ( '\\\\' . | ~ ( ( '\\\\' | '\\'' ) ) )* '\\'' )
            int alt11=2;
            int LA11_0 = input.LA(1);

            if ( (LA11_0=='\"') ) {
                alt11=1;
            }
            else if ( (LA11_0=='\'') ) {
                alt11=2;
            }
            else {
                NoViableAltException nvae =
                    new NoViableAltException("", 11, 0, input);

                throw nvae;
            }
            switch (alt11) {
                case 1 :
                    // InternalElkGraph.g:5995:16: '\"' ( '\\\\' . | ~ ( ( '\\\\' | '\"' ) ) )* '\"'
                    {
                    match('\"'); 
                    // InternalElkGraph.g:5995:20: ( '\\\\' . | ~ ( ( '\\\\' | '\"' ) ) )*
                    loop9:
                    do {
                        int alt9=3;
                        int LA9_0 = input.LA(1);

                        if ( (LA9_0=='\\') ) {
                            alt9=1;
                        }
                        else if ( ((LA9_0>='\u0000' && LA9_0<='!')||(LA9_0>='#' && LA9_0<='[')||(LA9_0>=']' && LA9_0<='\uFFFF')) ) {
                            alt9=2;
                        }


                        switch (alt9) {
                    	case 1 :
                    	    // InternalElkGraph.g:5995:21: '\\\\' .
                    	    {
                    	    match('\\'); 
                    	    matchAny(); 

                    	    }
                    	    break;
                    	case 2 :
                    	    // InternalElkGraph.g:5995:28: ~ ( ( '\\\\' | '\"' ) )
                    	    {
                    	    if ( (input.LA(1)>='\u0000' && input.LA(1)<='!')||(input.LA(1)>='#' && input.LA(1)<='[')||(input.LA(1)>=']' && input.LA(1)<='\uFFFF') ) {
                    	        input.consume();

                    	    }
                    	    else {
                    	        MismatchedSetException mse = new MismatchedSetException(null,input);
                    	        recover(mse);
                    	        throw mse;}


                    	    }
                    	    break;

                    	default :
                    	    break loop9;
                        }
                    } while (true);

                    match('\"'); 

                    }
                    break;
                case 2 :
                    // InternalElkGraph.g:5995:48: '\\'' ( '\\\\' . | ~ ( ( '\\\\' | '\\'' ) ) )* '\\''
                    {
                    match('\''); 
                    // InternalElkGraph.g:5995:53: ( '\\\\' . | ~ ( ( '\\\\' | '\\'' ) ) )*
                    loop10:
                    do {
                        int alt10=3;
                        int LA10_0 = input.LA(1);

                        if ( (LA10_0=='\\') ) {
                            alt10=1;
                        }
                        else if ( ((LA10_0>='\u0000' && LA10_0<='&')||(LA10_0>='(' && LA10_0<='[')||(LA10_0>=']' && LA10_0<='\uFFFF')) ) {
                            alt10=2;
                        }


                        switch (alt10) {
                    	case 1 :
                    	    // InternalElkGraph.g:5995:54: '\\\\' .
                    	    {
                    	    match('\\'); 
                    	    matchAny(); 

                    	    }
                    	    break;
                    	case 2 :
                    	    // InternalElkGraph.g:5995:61: ~ ( ( '\\\\' | '\\'' ) )
                    	    {
                    	    if ( (input.LA(1)>='\u0000' && input.LA(1)<='&')||(input.LA(1)>='(' && input.LA(1)<='[')||(input.LA(1)>=']' && input.LA(1)<='\uFFFF') ) {
                    	        input.consume();

                    	    }
                    	    else {
                    	        MismatchedSetException mse = new MismatchedSetException(null,input);
                    	        recover(mse);
                    	        throw mse;}


                    	    }
                    	    break;

                    	default :
                    	    break loop10;
                        }
                    } while (true);

                    match('\''); 

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
    // $ANTLR end "RULE_STRING"

    // $ANTLR start "RULE_ML_COMMENT"
    public final void mRULE_ML_COMMENT() throws RecognitionException {
        try {
            int _type = RULE_ML_COMMENT;
            int _channel = DEFAULT_TOKEN_CHANNEL;
            // InternalElkGraph.g:5997:17: ( '/*' ( options {greedy=false; } : . )* '*/' )
            // InternalElkGraph.g:5997:19: '/*' ( options {greedy=false; } : . )* '*/'
            {
            match("/*"); 

            // InternalElkGraph.g:5997:24: ( options {greedy=false; } : . )*
            loop12:
            do {
                int alt12=2;
                int LA12_0 = input.LA(1);

                if ( (LA12_0=='*') ) {
                    int LA12_1 = input.LA(2);

                    if ( (LA12_1=='/') ) {
                        alt12=2;
                    }
                    else if ( ((LA12_1>='\u0000' && LA12_1<='.')||(LA12_1>='0' && LA12_1<='\uFFFF')) ) {
                        alt12=1;
                    }


                }
                else if ( ((LA12_0>='\u0000' && LA12_0<=')')||(LA12_0>='+' && LA12_0<='\uFFFF')) ) {
                    alt12=1;
                }


                switch (alt12) {
            	case 1 :
            	    // InternalElkGraph.g:5997:52: .
            	    {
            	    matchAny(); 

            	    }
            	    break;

            	default :
            	    break loop12;
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
            // InternalElkGraph.g:5999:17: ( '//' (~ ( ( '\\n' | '\\r' ) ) )* ( ( '\\r' )? '\\n' )? )
            // InternalElkGraph.g:5999:19: '//' (~ ( ( '\\n' | '\\r' ) ) )* ( ( '\\r' )? '\\n' )?
            {
            match("//"); 

            // InternalElkGraph.g:5999:24: (~ ( ( '\\n' | '\\r' ) ) )*
            loop13:
            do {
                int alt13=2;
                int LA13_0 = input.LA(1);

                if ( ((LA13_0>='\u0000' && LA13_0<='\t')||(LA13_0>='\u000B' && LA13_0<='\f')||(LA13_0>='\u000E' && LA13_0<='\uFFFF')) ) {
                    alt13=1;
                }


                switch (alt13) {
            	case 1 :
            	    // InternalElkGraph.g:5999:24: ~ ( ( '\\n' | '\\r' ) )
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
            	    break loop13;
                }
            } while (true);

            // InternalElkGraph.g:5999:40: ( ( '\\r' )? '\\n' )?
            int alt15=2;
            int LA15_0 = input.LA(1);

            if ( (LA15_0=='\n'||LA15_0=='\r') ) {
                alt15=1;
            }
            switch (alt15) {
                case 1 :
                    // InternalElkGraph.g:5999:41: ( '\\r' )? '\\n'
                    {
                    // InternalElkGraph.g:5999:41: ( '\\r' )?
                    int alt14=2;
                    int LA14_0 = input.LA(1);

                    if ( (LA14_0=='\r') ) {
                        alt14=1;
                    }
                    switch (alt14) {
                        case 1 :
                            // InternalElkGraph.g:5999:41: '\\r'
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
            // InternalElkGraph.g:6001:9: ( ( ' ' | '\\t' | '\\r' | '\\n' )+ )
            // InternalElkGraph.g:6001:11: ( ' ' | '\\t' | '\\r' | '\\n' )+
            {
            // InternalElkGraph.g:6001:11: ( ' ' | '\\t' | '\\r' | '\\n' )+
            int cnt16=0;
            loop16:
            do {
                int alt16=2;
                int LA16_0 = input.LA(1);

                if ( ((LA16_0>='\t' && LA16_0<='\n')||LA16_0=='\r'||LA16_0==' ') ) {
                    alt16=1;
                }


                switch (alt16) {
            	case 1 :
            	    // InternalElkGraph.g:
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
            	    if ( cnt16 >= 1 ) break loop16;
                        EarlyExitException eee =
                            new EarlyExitException(16, input);
                        throw eee;
                }
                cnt16++;
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
            // InternalElkGraph.g:6003:16: ( . )
            // InternalElkGraph.g:6003:18: .
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
        // InternalElkGraph.g:1:8: ( T__13 | T__14 | T__15 | T__16 | T__17 | T__18 | T__19 | T__20 | T__21 | T__22 | T__23 | T__24 | T__25 | T__26 | T__27 | T__28 | T__29 | T__30 | T__31 | T__32 | T__33 | T__34 | T__35 | T__36 | T__37 | T__38 | RULE_SIGNED_INT | RULE_FLOAT | RULE_ID | RULE_STRING | RULE_ML_COMMENT | RULE_SL_COMMENT | RULE_WS | RULE_ANY_OTHER )
        int alt17=34;
        alt17 = dfa17.predict(input);
        switch (alt17) {
            case 1 :
                // InternalElkGraph.g:1:10: T__13
                {
                mT__13(); 

                }
                break;
            case 2 :
                // InternalElkGraph.g:1:16: T__14
                {
                mT__14(); 

                }
                break;
            case 3 :
                // InternalElkGraph.g:1:22: T__15
                {
                mT__15(); 

                }
                break;
            case 4 :
                // InternalElkGraph.g:1:28: T__16
                {
                mT__16(); 

                }
                break;
            case 5 :
                // InternalElkGraph.g:1:34: T__17
                {
                mT__17(); 

                }
                break;
            case 6 :
                // InternalElkGraph.g:1:40: T__18
                {
                mT__18(); 

                }
                break;
            case 7 :
                // InternalElkGraph.g:1:46: T__19
                {
                mT__19(); 

                }
                break;
            case 8 :
                // InternalElkGraph.g:1:52: T__20
                {
                mT__20(); 

                }
                break;
            case 9 :
                // InternalElkGraph.g:1:58: T__21
                {
                mT__21(); 

                }
                break;
            case 10 :
                // InternalElkGraph.g:1:64: T__22
                {
                mT__22(); 

                }
                break;
            case 11 :
                // InternalElkGraph.g:1:70: T__23
                {
                mT__23(); 

                }
                break;
            case 12 :
                // InternalElkGraph.g:1:76: T__24
                {
                mT__24(); 

                }
                break;
            case 13 :
                // InternalElkGraph.g:1:82: T__25
                {
                mT__25(); 

                }
                break;
            case 14 :
                // InternalElkGraph.g:1:88: T__26
                {
                mT__26(); 

                }
                break;
            case 15 :
                // InternalElkGraph.g:1:94: T__27
                {
                mT__27(); 

                }
                break;
            case 16 :
                // InternalElkGraph.g:1:100: T__28
                {
                mT__28(); 

                }
                break;
            case 17 :
                // InternalElkGraph.g:1:106: T__29
                {
                mT__29(); 

                }
                break;
            case 18 :
                // InternalElkGraph.g:1:112: T__30
                {
                mT__30(); 

                }
                break;
            case 19 :
                // InternalElkGraph.g:1:118: T__31
                {
                mT__31(); 

                }
                break;
            case 20 :
                // InternalElkGraph.g:1:124: T__32
                {
                mT__32(); 

                }
                break;
            case 21 :
                // InternalElkGraph.g:1:130: T__33
                {
                mT__33(); 

                }
                break;
            case 22 :
                // InternalElkGraph.g:1:136: T__34
                {
                mT__34(); 

                }
                break;
            case 23 :
                // InternalElkGraph.g:1:142: T__35
                {
                mT__35(); 

                }
                break;
            case 24 :
                // InternalElkGraph.g:1:148: T__36
                {
                mT__36(); 

                }
                break;
            case 25 :
                // InternalElkGraph.g:1:154: T__37
                {
                mT__37(); 

                }
                break;
            case 26 :
                // InternalElkGraph.g:1:160: T__38
                {
                mT__38(); 

                }
                break;
            case 27 :
                // InternalElkGraph.g:1:166: RULE_SIGNED_INT
                {
                mRULE_SIGNED_INT(); 

                }
                break;
            case 28 :
                // InternalElkGraph.g:1:182: RULE_FLOAT
                {
                mRULE_FLOAT(); 

                }
                break;
            case 29 :
                // InternalElkGraph.g:1:193: RULE_ID
                {
                mRULE_ID(); 

                }
                break;
            case 30 :
                // InternalElkGraph.g:1:201: RULE_STRING
                {
                mRULE_STRING(); 

                }
                break;
            case 31 :
                // InternalElkGraph.g:1:213: RULE_ML_COMMENT
                {
                mRULE_ML_COMMENT(); 

                }
                break;
            case 32 :
                // InternalElkGraph.g:1:229: RULE_SL_COMMENT
                {
                mRULE_SL_COMMENT(); 

                }
                break;
            case 33 :
                // InternalElkGraph.g:1:245: RULE_WS
                {
                mRULE_WS(); 

                }
                break;
            case 34 :
                // InternalElkGraph.g:1:253: RULE_ANY_OTHER
                {
                mRULE_ANY_OTHER(); 

                }
                break;

        }

    }


    protected DFA5 dfa5 = new DFA5(this);
    protected DFA17 dfa17 = new DFA17(this);
    static final String DFA5_eotS =
        "\4\uffff\1\5\1\uffff";
    static final String DFA5_eofS =
        "\6\uffff";
    static final String DFA5_minS =
        "\1\60\1\56\1\60\1\uffff\1\60\1\uffff";
    static final String DFA5_maxS =
        "\1\71\1\145\1\71\1\uffff\1\145\1\uffff";
    static final String DFA5_acceptS =
        "\3\uffff\1\2\1\uffff\1\1";
    static final String DFA5_specialS =
        "\6\uffff}>";
    static final String[] DFA5_transitionS = {
            "\12\1",
            "\1\2\1\uffff\12\1\13\uffff\1\3\37\uffff\1\3",
            "\12\4",
            "",
            "\12\4\13\uffff\1\3\37\uffff\1\3",
            ""
    };

    static final short[] DFA5_eot = DFA.unpackEncodedString(DFA5_eotS);
    static final short[] DFA5_eof = DFA.unpackEncodedString(DFA5_eofS);
    static final char[] DFA5_min = DFA.unpackEncodedStringToUnsignedChars(DFA5_minS);
    static final char[] DFA5_max = DFA.unpackEncodedStringToUnsignedChars(DFA5_maxS);
    static final short[] DFA5_accept = DFA.unpackEncodedString(DFA5_acceptS);
    static final short[] DFA5_special = DFA.unpackEncodedString(DFA5_specialS);
    static final short[][] DFA5_transition;

    static {
        int numStates = DFA5_transitionS.length;
        DFA5_transition = new short[numStates][];
        for (int i=0; i<numStates; i++) {
            DFA5_transition[i] = DFA.unpackEncodedString(DFA5_transitionS[i]);
        }
    }

    class DFA5 extends DFA {

        public DFA5(BaseRecognizer recognizer) {
            this.recognizer = recognizer;
            this.decisionNumber = 5;
            this.eot = DFA5_eot;
            this.eof = DFA5_eof;
            this.min = DFA5_min;
            this.max = DFA5_max;
            this.accept = DFA5_accept;
            this.special = DFA5_special;
            this.transition = DFA5_transition;
        }
        public String getDescription() {
            return "5989:25: ( RULE_INT '.' RULE_INT | RULE_INT ( '.' RULE_INT )? ( 'e' | 'E' ) ( '+' | '-' )? RULE_INT )";
        }
    }
    static final String DFA17_eotS =
        "\1\uffff\4\40\2\uffff\1\40\1\uffff\1\40\3\uffff\2\40\1\35\3\40\2\uffff\1\35\1\70\1\35\1\uffff\3\35\2\uffff\2\40\1\uffff\3\40\2\uffff\1\40\1\uffff\1\40\3\uffff\5\40\1\uffff\1\70\3\40\10\uffff\15\40\1\134\3\40\1\140\1\141\1\142\4\40\1\147\1\40\1\151\2\40\1\154\1\uffff\3\40\3\uffff\1\160\1\161\1\162\1\40\1\uffff\1\40\1\uffff\1\165\1\40\1\uffff\2\40\1\171\3\uffff\1\172\1\40\1\uffff\3\40\2\uffff\1\40\1\u0080\2\40\1\u0083\1\uffff\1\u0084\1\u0085\3\uffff";
    static final String DFA17_eofS =
        "\u0086\uffff";
    static final String DFA17_minS =
        "\1\0\1\157\1\162\1\141\1\162\2\uffff\1\141\1\uffff\1\157\3\uffff\1\145\1\144\1\60\1\156\1\165\1\145\2\uffff\1\60\1\56\1\101\1\uffff\2\0\1\52\2\uffff\1\154\1\144\1\uffff\1\165\1\154\1\141\2\uffff\1\142\1\uffff\1\162\3\uffff\1\172\1\141\1\143\1\147\1\144\1\uffff\1\56\1\143\1\164\1\156\10\uffff\1\154\2\145\1\163\1\160\1\145\1\157\1\164\1\151\1\145\1\162\1\164\1\145\1\60\1\157\1\147\1\144\3\60\1\145\1\150\1\154\1\165\1\60\1\164\1\60\1\164\1\151\1\60\1\uffff\1\155\1\157\1\163\3\uffff\3\60\1\164\1\uffff\1\151\1\uffff\1\60\1\157\1\uffff\2\151\1\60\3\uffff\1\60\1\157\1\uffff\3\156\2\uffff\1\156\1\60\2\147\1\60\1\uffff\2\60\3\uffff";
    static final String DFA17_maxS =
        "\1\uffff\1\165\1\162\1\141\1\162\2\uffff\1\141\1\uffff\1\157\3\uffff\1\164\1\156\1\76\1\156\1\165\1\145\2\uffff\1\71\1\145\1\172\1\uffff\2\uffff\1\57\2\uffff\1\154\1\144\1\uffff\1\165\1\154\1\141\2\uffff\1\171\1\uffff\1\163\3\uffff\1\172\1\141\1\143\1\147\1\144\1\uffff\1\145\1\143\1\164\1\156\10\uffff\1\154\2\145\1\163\1\160\1\145\1\157\1\164\1\151\1\145\1\162\1\164\1\145\1\172\1\157\1\147\1\144\3\172\1\145\1\150\1\154\1\165\1\172\1\164\1\172\1\164\1\151\1\172\1\uffff\1\155\1\157\1\163\3\uffff\3\172\1\164\1\uffff\1\151\1\uffff\1\172\1\157\1\uffff\2\151\1\172\3\uffff\1\172\1\157\1\uffff\3\156\2\uffff\1\156\1\172\2\147\1\172\1\uffff\2\172\3\uffff";
    static final String DFA17_acceptS =
        "\5\uffff\1\6\1\7\1\uffff\1\11\1\uffff\1\14\1\15\1\17\6\uffff\1\30\1\32\3\uffff\1\35\3\uffff\1\41\1\42\2\uffff\1\35\3\uffff\1\6\1\7\1\uffff\1\11\1\uffff\1\14\1\15\1\17\5\uffff\1\22\4\uffff\1\30\1\32\1\33\1\34\1\36\1\37\1\40\1\41\36\uffff\1\26\3\uffff\1\1\1\5\1\2\4\uffff\1\12\1\uffff\1\20\2\uffff\1\21\3\uffff\1\3\1\4\1\10\2\uffff\1\25\3\uffff\1\27\1\13\5\uffff\1\31\2\uffff\1\16\1\23\1\24";
    static final String DFA17_specialS =
        "\1\0\30\uffff\1\1\1\2\153\uffff}>";
    static final String[] DFA17_transitionS = {
            "\11\35\2\34\2\35\1\34\22\35\1\34\1\35\1\31\4\35\1\32\3\35\1\25\1\14\1\17\1\24\1\33\12\26\1\10\6\35\32\30\1\12\1\35\1\13\1\27\1\30\1\35\1\30\1\22\2\30\1\16\1\3\1\4\1\30\1\20\2\30\1\7\1\30\1\1\1\21\1\11\2\30\1\15\1\2\6\30\1\5\1\23\1\6\uff82\35",
            "\1\37\5\uffff\1\36",
            "\1\41",
            "\1\42",
            "\1\43",
            "",
            "",
            "\1\46",
            "",
            "\1\50",
            "",
            "",
            "",
            "\1\56\3\uffff\1\54\12\uffff\1\55",
            "\1\57\11\uffff\1\60",
            "\12\62\4\uffff\1\61",
            "\1\63",
            "\1\64",
            "\1\65",
            "",
            "",
            "\12\62",
            "\1\71\1\uffff\12\62\13\uffff\1\71\37\uffff\1\71",
            "\32\40\4\uffff\1\40\1\uffff\32\40",
            "",
            "\0\72",
            "\0\72",
            "\1\73\4\uffff\1\74",
            "",
            "",
            "\1\76",
            "\1\77",
            "",
            "\1\100",
            "\1\101",
            "\1\102",
            "",
            "",
            "\1\103\26\uffff\1\104",
            "",
            "\1\105\1\106",
            "",
            "",
            "",
            "\1\107",
            "\1\110",
            "\1\111",
            "\1\112",
            "\1\113",
            "",
            "\1\71\1\uffff\12\62\13\uffff\1\71\37\uffff\1\71",
            "\1\114",
            "\1\115",
            "\1\116",
            "",
            "",
            "",
            "",
            "",
            "",
            "",
            "",
            "\1\117",
            "\1\120",
            "\1\121",
            "\1\122",
            "\1\123",
            "\1\124",
            "\1\125",
            "\1\126",
            "\1\127",
            "\1\130",
            "\1\131",
            "\1\132",
            "\1\133",
            "\12\40\7\uffff\32\40\4\uffff\1\40\1\uffff\32\40",
            "\1\135",
            "\1\136",
            "\1\137",
            "\12\40\7\uffff\32\40\4\uffff\1\40\1\uffff\32\40",
            "\12\40\7\uffff\32\40\4\uffff\1\40\1\uffff\32\40",
            "\12\40\7\uffff\32\40\4\uffff\1\40\1\uffff\32\40",
            "\1\143",
            "\1\144",
            "\1\145",
            "\1\146",
            "\12\40\7\uffff\32\40\4\uffff\1\40\1\uffff\32\40",
            "\1\150",
            "\12\40\7\uffff\32\40\4\uffff\1\40\1\uffff\32\40",
            "\1\152",
            "\1\153",
            "\12\40\7\uffff\32\40\4\uffff\1\40\1\uffff\32\40",
            "",
            "\1\155",
            "\1\156",
            "\1\157",
            "",
            "",
            "",
            "\12\40\7\uffff\32\40\4\uffff\1\40\1\uffff\32\40",
            "\12\40\7\uffff\32\40\4\uffff\1\40\1\uffff\32\40",
            "\12\40\7\uffff\32\40\4\uffff\1\40\1\uffff\32\40",
            "\1\163",
            "",
            "\1\164",
            "",
            "\12\40\7\uffff\32\40\4\uffff\1\40\1\uffff\32\40",
            "\1\166",
            "",
            "\1\167",
            "\1\170",
            "\12\40\7\uffff\32\40\4\uffff\1\40\1\uffff\32\40",
            "",
            "",
            "",
            "\12\40\7\uffff\32\40\4\uffff\1\40\1\uffff\32\40",
            "\1\173",
            "",
            "\1\174",
            "\1\175",
            "\1\176",
            "",
            "",
            "\1\177",
            "\12\40\7\uffff\32\40\4\uffff\1\40\1\uffff\32\40",
            "\1\u0081",
            "\1\u0082",
            "\12\40\7\uffff\32\40\4\uffff\1\40\1\uffff\32\40",
            "",
            "\12\40\7\uffff\32\40\4\uffff\1\40\1\uffff\32\40",
            "\12\40\7\uffff\32\40\4\uffff\1\40\1\uffff\32\40",
            "",
            "",
            ""
    };

    static final short[] DFA17_eot = DFA.unpackEncodedString(DFA17_eotS);
    static final short[] DFA17_eof = DFA.unpackEncodedString(DFA17_eofS);
    static final char[] DFA17_min = DFA.unpackEncodedStringToUnsignedChars(DFA17_minS);
    static final char[] DFA17_max = DFA.unpackEncodedStringToUnsignedChars(DFA17_maxS);
    static final short[] DFA17_accept = DFA.unpackEncodedString(DFA17_acceptS);
    static final short[] DFA17_special = DFA.unpackEncodedString(DFA17_specialS);
    static final short[][] DFA17_transition;

    static {
        int numStates = DFA17_transitionS.length;
        DFA17_transition = new short[numStates][];
        for (int i=0; i<numStates; i++) {
            DFA17_transition[i] = DFA.unpackEncodedString(DFA17_transitionS[i]);
        }
    }

    class DFA17 extends DFA {

        public DFA17(BaseRecognizer recognizer) {
            this.recognizer = recognizer;
            this.decisionNumber = 17;
            this.eot = DFA17_eot;
            this.eof = DFA17_eof;
            this.min = DFA17_min;
            this.max = DFA17_max;
            this.accept = DFA17_accept;
            this.special = DFA17_special;
            this.transition = DFA17_transition;
        }
        public String getDescription() {
            return "1:1: Tokens : ( T__13 | T__14 | T__15 | T__16 | T__17 | T__18 | T__19 | T__20 | T__21 | T__22 | T__23 | T__24 | T__25 | T__26 | T__27 | T__28 | T__29 | T__30 | T__31 | T__32 | T__33 | T__34 | T__35 | T__36 | T__37 | T__38 | RULE_SIGNED_INT | RULE_FLOAT | RULE_ID | RULE_STRING | RULE_ML_COMMENT | RULE_SL_COMMENT | RULE_WS | RULE_ANY_OTHER );";
        }
        public int specialStateTransition(int s, IntStream _input) throws NoViableAltException {
            IntStream input = _input;
        	int _s = s;
            switch ( s ) {
                    case 0 : 
                        int LA17_0 = input.LA(1);

                        s = -1;
                        if ( (LA17_0=='n') ) {s = 1;}

                        else if ( (LA17_0=='t') ) {s = 2;}

                        else if ( (LA17_0=='f') ) {s = 3;}

                        else if ( (LA17_0=='g') ) {s = 4;}

                        else if ( (LA17_0=='{') ) {s = 5;}

                        else if ( (LA17_0=='}') ) {s = 6;}

                        else if ( (LA17_0=='l') ) {s = 7;}

                        else if ( (LA17_0==':') ) {s = 8;}

                        else if ( (LA17_0=='p') ) {s = 9;}

                        else if ( (LA17_0=='[') ) {s = 10;}

                        else if ( (LA17_0==']') ) {s = 11;}

                        else if ( (LA17_0==',') ) {s = 12;}

                        else if ( (LA17_0=='s') ) {s = 13;}

                        else if ( (LA17_0=='e') ) {s = 14;}

                        else if ( (LA17_0=='-') ) {s = 15;}

                        else if ( (LA17_0=='i') ) {s = 16;}

                        else if ( (LA17_0=='o') ) {s = 17;}

                        else if ( (LA17_0=='b') ) {s = 18;}

                        else if ( (LA17_0=='|') ) {s = 19;}

                        else if ( (LA17_0=='.') ) {s = 20;}

                        else if ( (LA17_0=='+') ) {s = 21;}

                        else if ( ((LA17_0>='0' && LA17_0<='9')) ) {s = 22;}

                        else if ( (LA17_0=='^') ) {s = 23;}

                        else if ( ((LA17_0>='A' && LA17_0<='Z')||LA17_0=='_'||LA17_0=='a'||(LA17_0>='c' && LA17_0<='d')||LA17_0=='h'||(LA17_0>='j' && LA17_0<='k')||LA17_0=='m'||(LA17_0>='q' && LA17_0<='r')||(LA17_0>='u' && LA17_0<='z')) ) {s = 24;}

                        else if ( (LA17_0=='\"') ) {s = 25;}

                        else if ( (LA17_0=='\'') ) {s = 26;}

                        else if ( (LA17_0=='/') ) {s = 27;}

                        else if ( ((LA17_0>='\t' && LA17_0<='\n')||LA17_0=='\r'||LA17_0==' ') ) {s = 28;}

                        else if ( ((LA17_0>='\u0000' && LA17_0<='\b')||(LA17_0>='\u000B' && LA17_0<='\f')||(LA17_0>='\u000E' && LA17_0<='\u001F')||LA17_0=='!'||(LA17_0>='#' && LA17_0<='&')||(LA17_0>='(' && LA17_0<='*')||(LA17_0>=';' && LA17_0<='@')||LA17_0=='\\'||LA17_0=='`'||(LA17_0>='~' && LA17_0<='\uFFFF')) ) {s = 29;}

                        if ( s>=0 ) return s;
                        break;
                    case 1 : 
                        int LA17_25 = input.LA(1);

                        s = -1;
                        if ( ((LA17_25>='\u0000' && LA17_25<='\uFFFF')) ) {s = 58;}

                        else s = 29;

                        if ( s>=0 ) return s;
                        break;
                    case 2 : 
                        int LA17_26 = input.LA(1);

                        s = -1;
                        if ( ((LA17_26>='\u0000' && LA17_26<='\uFFFF')) ) {s = 58;}

                        else s = 29;

                        if ( s>=0 ) return s;
                        break;
            }
            NoViableAltException nvae =
                new NoViableAltException(getDescription(), 17, _s, input);
            error(nvae);
            throw nvae;
        }
    }
 

}