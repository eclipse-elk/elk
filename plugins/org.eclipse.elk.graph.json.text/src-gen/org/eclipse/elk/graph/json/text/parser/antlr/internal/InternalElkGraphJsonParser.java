package org.eclipse.elk.graph.json.text.parser.antlr.internal;

import org.eclipse.xtext.*;
import org.eclipse.xtext.parser.*;
import org.eclipse.xtext.parser.impl.*;
import org.eclipse.emf.ecore.util.EcoreUtil;
import org.eclipse.emf.ecore.EObject;
import org.eclipse.xtext.parser.antlr.AbstractInternalAntlrParser;
import org.eclipse.xtext.parser.antlr.XtextTokenStream;
import org.eclipse.xtext.parser.antlr.XtextTokenStream.HiddenTokens;
import org.eclipse.xtext.parser.antlr.AntlrDatatypeRuleToken;
import org.eclipse.elk.graph.json.text.services.ElkGraphJsonGrammarAccess;



import org.antlr.runtime.*;
import java.util.Stack;
import java.util.List;
import java.util.ArrayList;

/*******************************************************************************
 * Copyright (c) 2020 Kiel University and others.
 *
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
@SuppressWarnings("all")
public class InternalElkGraphJsonParser extends AbstractInternalAntlrParser {
    public static final String[] tokenNames = new String[] {
        "<invalid>", "<EOR>", "<DOWN>", "<UP>", "RULE_STRING", "RULE_ID", "RULE_SIGNED_INT", "RULE_FLOAT", "RULE_INT", "RULE_ML_COMMENT", "RULE_SL_COMMENT", "RULE_WS", "RULE_ANY_OTHER", "'{'", "','", "'}'", "':'", "'['", "']'", "'true'", "'false'", "'null'", "'\"children\"'", "'\\'children\\''", "'children'", "'\"ports\"'", "'\\'ports\\''", "'ports'", "'\"labels\"'", "'\\'labels\\''", "'labels'", "'\"edges\"'", "'\\'edges\\''", "'edges'", "'\"layoutOptions\"'", "'\\'layoutOptions\\''", "'layoutOptions'", "'\"properties\"'", "'\\'properties\\''", "'properties'", "'\"id\"'", "'\\'id\\''", "'id'", "'\"x\"'", "'\\'x\\''", "'x'", "'\"y\"'", "'\\'y\\''", "'y'", "'\"width\"'", "'\\'width\\''", "'width'", "'\"height\"'", "'\\'height\\''", "'height'", "'\"sources\"'", "'\\'sources\\''", "'sources'", "'\"targets\"'", "'\\'targets\\''", "'targets'", "'\"text\"'", "'\\'text\\''", "'text'"
    };
    public static final int T__50=50;
    public static final int T__19=19;
    public static final int T__15=15;
    public static final int T__59=59;
    public static final int T__16=16;
    public static final int T__17=17;
    public static final int T__18=18;
    public static final int T__55=55;
    public static final int T__56=56;
    public static final int T__13=13;
    public static final int T__57=57;
    public static final int T__14=14;
    public static final int T__58=58;
    public static final int T__51=51;
    public static final int T__52=52;
    public static final int T__53=53;
    public static final int T__54=54;
    public static final int T__60=60;
    public static final int T__61=61;
    public static final int RULE_ID=5;
    public static final int T__26=26;
    public static final int T__27=27;
    public static final int T__28=28;
    public static final int RULE_INT=8;
    public static final int T__29=29;
    public static final int T__22=22;
    public static final int RULE_ML_COMMENT=9;
    public static final int T__23=23;
    public static final int T__24=24;
    public static final int T__25=25;
    public static final int T__62=62;
    public static final int T__63=63;
    public static final int T__20=20;
    public static final int T__21=21;
    public static final int RULE_STRING=4;
    public static final int RULE_SL_COMMENT=10;
    public static final int T__37=37;
    public static final int T__38=38;
    public static final int T__39=39;
    public static final int T__33=33;
    public static final int T__34=34;
    public static final int T__35=35;
    public static final int T__36=36;
    public static final int EOF=-1;
    public static final int T__30=30;
    public static final int T__31=31;
    public static final int T__32=32;
    public static final int RULE_WS=11;
    public static final int RULE_SIGNED_INT=6;
    public static final int RULE_ANY_OTHER=12;
    public static final int T__48=48;
    public static final int T__49=49;
    public static final int T__44=44;
    public static final int T__45=45;
    public static final int RULE_FLOAT=7;
    public static final int T__46=46;
    public static final int T__47=47;
    public static final int T__40=40;
    public static final int T__41=41;
    public static final int T__42=42;
    public static final int T__43=43;

    // delegates
    // delegators


        public InternalElkGraphJsonParser(TokenStream input) {
            this(input, new RecognizerSharedState());
        }
        public InternalElkGraphJsonParser(TokenStream input, RecognizerSharedState state) {
            super(input, state);
             
        }
        

    public String[] getTokenNames() { return InternalElkGraphJsonParser.tokenNames; }
    public String getGrammarFileName() { return "InternalElkGraphJson.g"; }



     	private ElkGraphJsonGrammarAccess grammarAccess;

        public InternalElkGraphJsonParser(TokenStream input, ElkGraphJsonGrammarAccess grammarAccess) {
            this(input);
            this.grammarAccess = grammarAccess;
            registerRules(grammarAccess.getGrammar());
        }

        @Override
        protected String getFirstRuleName() {
        	return "ElkNode";
       	}

       	@Override
       	protected ElkGraphJsonGrammarAccess getGrammarAccess() {
       		return grammarAccess;
       	}




    // $ANTLR start "entryRuleElkNode"
    // InternalElkGraphJson.g:70:1: entryRuleElkNode returns [EObject current=null] : iv_ruleElkNode= ruleElkNode EOF ;
    public final EObject entryRuleElkNode() throws RecognitionException {
        EObject current = null;

        EObject iv_ruleElkNode = null;


        try {
            // InternalElkGraphJson.g:70:48: (iv_ruleElkNode= ruleElkNode EOF )
            // InternalElkGraphJson.g:71:2: iv_ruleElkNode= ruleElkNode EOF
            {
             newCompositeNode(grammarAccess.getElkNodeRule()); 
            pushFollow(FOLLOW_1);
            iv_ruleElkNode=ruleElkNode();

            state._fsp--;

             current =iv_ruleElkNode; 
            match(input,EOF,FOLLOW_2); 

            }

        }

            catch (RecognitionException re) {
                recover(input,re);
                appendSkippedTokens();
            }
        finally {
        }
        return current;
    }
    // $ANTLR end "entryRuleElkNode"


    // $ANTLR start "ruleElkNode"
    // InternalElkGraphJson.g:77:1: ruleElkNode returns [EObject current=null] : ( () otherlv_1= '{' (this_NodeElement_2= ruleNodeElement[$current] (otherlv_3= ',' this_NodeElement_4= ruleNodeElement[$current] )* )? (otherlv_5= ',' )? otherlv_6= '}' ) ;
    public final EObject ruleElkNode() throws RecognitionException {
        EObject current = null;

        Token otherlv_1=null;
        Token otherlv_3=null;
        Token otherlv_5=null;
        Token otherlv_6=null;
        EObject this_NodeElement_2 = null;

        EObject this_NodeElement_4 = null;



        	enterRule();

        try {
            // InternalElkGraphJson.g:83:2: ( ( () otherlv_1= '{' (this_NodeElement_2= ruleNodeElement[$current] (otherlv_3= ',' this_NodeElement_4= ruleNodeElement[$current] )* )? (otherlv_5= ',' )? otherlv_6= '}' ) )
            // InternalElkGraphJson.g:84:2: ( () otherlv_1= '{' (this_NodeElement_2= ruleNodeElement[$current] (otherlv_3= ',' this_NodeElement_4= ruleNodeElement[$current] )* )? (otherlv_5= ',' )? otherlv_6= '}' )
            {
            // InternalElkGraphJson.g:84:2: ( () otherlv_1= '{' (this_NodeElement_2= ruleNodeElement[$current] (otherlv_3= ',' this_NodeElement_4= ruleNodeElement[$current] )* )? (otherlv_5= ',' )? otherlv_6= '}' )
            // InternalElkGraphJson.g:85:3: () otherlv_1= '{' (this_NodeElement_2= ruleNodeElement[$current] (otherlv_3= ',' this_NodeElement_4= ruleNodeElement[$current] )* )? (otherlv_5= ',' )? otherlv_6= '}'
            {
            // InternalElkGraphJson.g:85:3: ()
            // InternalElkGraphJson.g:86:4: 
            {

            				current = forceCreateModelElement(
            					grammarAccess.getElkNodeAccess().getElkNodeAction_0(),
            					current);
            			

            }

            otherlv_1=(Token)match(input,13,FOLLOW_3); 

            			newLeafNode(otherlv_1, grammarAccess.getElkNodeAccess().getLeftCurlyBracketKeyword_1());
            		
            // InternalElkGraphJson.g:96:3: (this_NodeElement_2= ruleNodeElement[$current] (otherlv_3= ',' this_NodeElement_4= ruleNodeElement[$current] )* )?
            int alt2=2;
            int LA2_0 = input.LA(1);

            if ( ((LA2_0>=RULE_STRING && LA2_0<=RULE_ID)||(LA2_0>=22 && LA2_0<=54)) ) {
                alt2=1;
            }
            switch (alt2) {
                case 1 :
                    // InternalElkGraphJson.g:97:4: this_NodeElement_2= ruleNodeElement[$current] (otherlv_3= ',' this_NodeElement_4= ruleNodeElement[$current] )*
                    {

                    				if (current==null) {
                    					current = createModelElement(grammarAccess.getElkNodeRule());
                    				}
                    				newCompositeNode(grammarAccess.getElkNodeAccess().getNodeElementParserRuleCall_2_0());
                    			
                    pushFollow(FOLLOW_4);
                    this_NodeElement_2=ruleNodeElement(current);

                    state._fsp--;


                    				current = this_NodeElement_2;
                    				afterParserOrEnumRuleCall();
                    			
                    // InternalElkGraphJson.g:108:4: (otherlv_3= ',' this_NodeElement_4= ruleNodeElement[$current] )*
                    loop1:
                    do {
                        int alt1=2;
                        int LA1_0 = input.LA(1);

                        if ( (LA1_0==14) ) {
                            int LA1_1 = input.LA(2);

                            if ( ((LA1_1>=RULE_STRING && LA1_1<=RULE_ID)||(LA1_1>=22 && LA1_1<=54)) ) {
                                alt1=1;
                            }


                        }


                        switch (alt1) {
                    	case 1 :
                    	    // InternalElkGraphJson.g:109:5: otherlv_3= ',' this_NodeElement_4= ruleNodeElement[$current]
                    	    {
                    	    otherlv_3=(Token)match(input,14,FOLLOW_5); 

                    	    					newLeafNode(otherlv_3, grammarAccess.getElkNodeAccess().getCommaKeyword_2_1_0());
                    	    				

                    	    					if (current==null) {
                    	    						current = createModelElement(grammarAccess.getElkNodeRule());
                    	    					}
                    	    					newCompositeNode(grammarAccess.getElkNodeAccess().getNodeElementParserRuleCall_2_1_1());
                    	    				
                    	    pushFollow(FOLLOW_4);
                    	    this_NodeElement_4=ruleNodeElement(current);

                    	    state._fsp--;


                    	    					current = this_NodeElement_4;
                    	    					afterParserOrEnumRuleCall();
                    	    				

                    	    }
                    	    break;

                    	default :
                    	    break loop1;
                        }
                    } while (true);


                    }
                    break;

            }

            // InternalElkGraphJson.g:126:3: (otherlv_5= ',' )?
            int alt3=2;
            int LA3_0 = input.LA(1);

            if ( (LA3_0==14) ) {
                alt3=1;
            }
            switch (alt3) {
                case 1 :
                    // InternalElkGraphJson.g:127:4: otherlv_5= ','
                    {
                    otherlv_5=(Token)match(input,14,FOLLOW_6); 

                    				newLeafNode(otherlv_5, grammarAccess.getElkNodeAccess().getCommaKeyword_3());
                    			

                    }
                    break;

            }

            otherlv_6=(Token)match(input,15,FOLLOW_2); 

            			newLeafNode(otherlv_6, grammarAccess.getElkNodeAccess().getRightCurlyBracketKeyword_4());
            		

            }


            }


            	leaveRule();

        }

            catch (RecognitionException re) {
                recover(input,re);
                appendSkippedTokens();
            }
        finally {
        }
        return current;
    }
    // $ANTLR end "ruleElkNode"


    // $ANTLR start "ruleNodeElement"
    // InternalElkGraphJson.g:141:1: ruleNodeElement[EObject in_current] returns [EObject current=in_current] : (this_ElkId_0= ruleElkId[$current] | this_ShapeElement_1= ruleShapeElement[$current] | ( ruleKeyChildren otherlv_3= ':' this_ElkNodeChildren_4= ruleElkNodeChildren[$current] ) | ( ruleKeyPorts otherlv_6= ':' this_ElkNodePorts_7= ruleElkNodePorts[$current] ) | ( ruleKeyLabels otherlv_9= ':' this_ElkGraphElementLabels_10= ruleElkGraphElementLabels[$current] ) | ( ruleKeyEdges otherlv_12= ':' this_ElkNodeEdges_13= ruleElkNodeEdges[$current] ) | ( ruleKeyLayoutOptions otherlv_15= ':' this_ElkGraphElementProperties_16= ruleElkGraphElementProperties[$current] ) | ruleJsonMember ) ;
    public final EObject ruleNodeElement(EObject in_current) throws RecognitionException {
        EObject current = in_current;

        Token otherlv_3=null;
        Token otherlv_6=null;
        Token otherlv_9=null;
        Token otherlv_12=null;
        Token otherlv_15=null;
        EObject this_ElkId_0 = null;

        EObject this_ShapeElement_1 = null;

        EObject this_ElkNodeChildren_4 = null;

        EObject this_ElkNodePorts_7 = null;

        EObject this_ElkGraphElementLabels_10 = null;

        EObject this_ElkNodeEdges_13 = null;

        EObject this_ElkGraphElementProperties_16 = null;



        	enterRule();

        try {
            // InternalElkGraphJson.g:147:2: ( (this_ElkId_0= ruleElkId[$current] | this_ShapeElement_1= ruleShapeElement[$current] | ( ruleKeyChildren otherlv_3= ':' this_ElkNodeChildren_4= ruleElkNodeChildren[$current] ) | ( ruleKeyPorts otherlv_6= ':' this_ElkNodePorts_7= ruleElkNodePorts[$current] ) | ( ruleKeyLabels otherlv_9= ':' this_ElkGraphElementLabels_10= ruleElkGraphElementLabels[$current] ) | ( ruleKeyEdges otherlv_12= ':' this_ElkNodeEdges_13= ruleElkNodeEdges[$current] ) | ( ruleKeyLayoutOptions otherlv_15= ':' this_ElkGraphElementProperties_16= ruleElkGraphElementProperties[$current] ) | ruleJsonMember ) )
            // InternalElkGraphJson.g:148:2: (this_ElkId_0= ruleElkId[$current] | this_ShapeElement_1= ruleShapeElement[$current] | ( ruleKeyChildren otherlv_3= ':' this_ElkNodeChildren_4= ruleElkNodeChildren[$current] ) | ( ruleKeyPorts otherlv_6= ':' this_ElkNodePorts_7= ruleElkNodePorts[$current] ) | ( ruleKeyLabels otherlv_9= ':' this_ElkGraphElementLabels_10= ruleElkGraphElementLabels[$current] ) | ( ruleKeyEdges otherlv_12= ':' this_ElkNodeEdges_13= ruleElkNodeEdges[$current] ) | ( ruleKeyLayoutOptions otherlv_15= ':' this_ElkGraphElementProperties_16= ruleElkGraphElementProperties[$current] ) | ruleJsonMember )
            {
            // InternalElkGraphJson.g:148:2: (this_ElkId_0= ruleElkId[$current] | this_ShapeElement_1= ruleShapeElement[$current] | ( ruleKeyChildren otherlv_3= ':' this_ElkNodeChildren_4= ruleElkNodeChildren[$current] ) | ( ruleKeyPorts otherlv_6= ':' this_ElkNodePorts_7= ruleElkNodePorts[$current] ) | ( ruleKeyLabels otherlv_9= ':' this_ElkGraphElementLabels_10= ruleElkGraphElementLabels[$current] ) | ( ruleKeyEdges otherlv_12= ':' this_ElkNodeEdges_13= ruleElkNodeEdges[$current] ) | ( ruleKeyLayoutOptions otherlv_15= ':' this_ElkGraphElementProperties_16= ruleElkGraphElementProperties[$current] ) | ruleJsonMember )
            int alt4=8;
            switch ( input.LA(1) ) {
            case 40:
            case 41:
            case 42:
                {
                alt4=1;
                }
                break;
            case 43:
            case 44:
            case 45:
            case 46:
            case 47:
            case 48:
            case 49:
            case 50:
            case 51:
            case 52:
            case 53:
            case 54:
                {
                alt4=2;
                }
                break;
            case 22:
            case 23:
            case 24:
                {
                alt4=3;
                }
                break;
            case 25:
            case 26:
            case 27:
                {
                alt4=4;
                }
                break;
            case 28:
            case 29:
            case 30:
                {
                alt4=5;
                }
                break;
            case 31:
            case 32:
            case 33:
                {
                alt4=6;
                }
                break;
            case 34:
            case 35:
            case 36:
            case 37:
            case 38:
            case 39:
                {
                alt4=7;
                }
                break;
            case RULE_STRING:
            case RULE_ID:
                {
                alt4=8;
                }
                break;
            default:
                NoViableAltException nvae =
                    new NoViableAltException("", 4, 0, input);

                throw nvae;
            }

            switch (alt4) {
                case 1 :
                    // InternalElkGraphJson.g:149:3: this_ElkId_0= ruleElkId[$current]
                    {

                    			if (current==null) {
                    				current = createModelElement(grammarAccess.getNodeElementRule());
                    			}
                    			newCompositeNode(grammarAccess.getNodeElementAccess().getElkIdParserRuleCall_0());
                    		
                    pushFollow(FOLLOW_2);
                    this_ElkId_0=ruleElkId(current);

                    state._fsp--;


                    			current = this_ElkId_0;
                    			afterParserOrEnumRuleCall();
                    		

                    }
                    break;
                case 2 :
                    // InternalElkGraphJson.g:161:3: this_ShapeElement_1= ruleShapeElement[$current]
                    {

                    			if (current==null) {
                    				current = createModelElement(grammarAccess.getNodeElementRule());
                    			}
                    			newCompositeNode(grammarAccess.getNodeElementAccess().getShapeElementParserRuleCall_1());
                    		
                    pushFollow(FOLLOW_2);
                    this_ShapeElement_1=ruleShapeElement(current);

                    state._fsp--;


                    			current = this_ShapeElement_1;
                    			afterParserOrEnumRuleCall();
                    		

                    }
                    break;
                case 3 :
                    // InternalElkGraphJson.g:173:3: ( ruleKeyChildren otherlv_3= ':' this_ElkNodeChildren_4= ruleElkNodeChildren[$current] )
                    {
                    // InternalElkGraphJson.g:173:3: ( ruleKeyChildren otherlv_3= ':' this_ElkNodeChildren_4= ruleElkNodeChildren[$current] )
                    // InternalElkGraphJson.g:174:4: ruleKeyChildren otherlv_3= ':' this_ElkNodeChildren_4= ruleElkNodeChildren[$current]
                    {

                    				newCompositeNode(grammarAccess.getNodeElementAccess().getKeyChildrenParserRuleCall_2_0());
                    			
                    pushFollow(FOLLOW_7);
                    ruleKeyChildren();

                    state._fsp--;


                    				afterParserOrEnumRuleCall();
                    			
                    otherlv_3=(Token)match(input,16,FOLLOW_8); 

                    				newLeafNode(otherlv_3, grammarAccess.getNodeElementAccess().getColonKeyword_2_1());
                    			

                    				if (current==null) {
                    					current = createModelElement(grammarAccess.getNodeElementRule());
                    				}
                    				newCompositeNode(grammarAccess.getNodeElementAccess().getElkNodeChildrenParserRuleCall_2_2());
                    			
                    pushFollow(FOLLOW_2);
                    this_ElkNodeChildren_4=ruleElkNodeChildren(current);

                    state._fsp--;


                    				current = this_ElkNodeChildren_4;
                    				afterParserOrEnumRuleCall();
                    			

                    }


                    }
                    break;
                case 4 :
                    // InternalElkGraphJson.g:198:3: ( ruleKeyPorts otherlv_6= ':' this_ElkNodePorts_7= ruleElkNodePorts[$current] )
                    {
                    // InternalElkGraphJson.g:198:3: ( ruleKeyPorts otherlv_6= ':' this_ElkNodePorts_7= ruleElkNodePorts[$current] )
                    // InternalElkGraphJson.g:199:4: ruleKeyPorts otherlv_6= ':' this_ElkNodePorts_7= ruleElkNodePorts[$current]
                    {

                    				newCompositeNode(grammarAccess.getNodeElementAccess().getKeyPortsParserRuleCall_3_0());
                    			
                    pushFollow(FOLLOW_7);
                    ruleKeyPorts();

                    state._fsp--;


                    				afterParserOrEnumRuleCall();
                    			
                    otherlv_6=(Token)match(input,16,FOLLOW_8); 

                    				newLeafNode(otherlv_6, grammarAccess.getNodeElementAccess().getColonKeyword_3_1());
                    			

                    				if (current==null) {
                    					current = createModelElement(grammarAccess.getNodeElementRule());
                    				}
                    				newCompositeNode(grammarAccess.getNodeElementAccess().getElkNodePortsParserRuleCall_3_2());
                    			
                    pushFollow(FOLLOW_2);
                    this_ElkNodePorts_7=ruleElkNodePorts(current);

                    state._fsp--;


                    				current = this_ElkNodePorts_7;
                    				afterParserOrEnumRuleCall();
                    			

                    }


                    }
                    break;
                case 5 :
                    // InternalElkGraphJson.g:223:3: ( ruleKeyLabels otherlv_9= ':' this_ElkGraphElementLabels_10= ruleElkGraphElementLabels[$current] )
                    {
                    // InternalElkGraphJson.g:223:3: ( ruleKeyLabels otherlv_9= ':' this_ElkGraphElementLabels_10= ruleElkGraphElementLabels[$current] )
                    // InternalElkGraphJson.g:224:4: ruleKeyLabels otherlv_9= ':' this_ElkGraphElementLabels_10= ruleElkGraphElementLabels[$current]
                    {

                    				newCompositeNode(grammarAccess.getNodeElementAccess().getKeyLabelsParserRuleCall_4_0());
                    			
                    pushFollow(FOLLOW_7);
                    ruleKeyLabels();

                    state._fsp--;


                    				afterParserOrEnumRuleCall();
                    			
                    otherlv_9=(Token)match(input,16,FOLLOW_8); 

                    				newLeafNode(otherlv_9, grammarAccess.getNodeElementAccess().getColonKeyword_4_1());
                    			

                    				if (current==null) {
                    					current = createModelElement(grammarAccess.getNodeElementRule());
                    				}
                    				newCompositeNode(grammarAccess.getNodeElementAccess().getElkGraphElementLabelsParserRuleCall_4_2());
                    			
                    pushFollow(FOLLOW_2);
                    this_ElkGraphElementLabels_10=ruleElkGraphElementLabels(current);

                    state._fsp--;


                    				current = this_ElkGraphElementLabels_10;
                    				afterParserOrEnumRuleCall();
                    			

                    }


                    }
                    break;
                case 6 :
                    // InternalElkGraphJson.g:248:3: ( ruleKeyEdges otherlv_12= ':' this_ElkNodeEdges_13= ruleElkNodeEdges[$current] )
                    {
                    // InternalElkGraphJson.g:248:3: ( ruleKeyEdges otherlv_12= ':' this_ElkNodeEdges_13= ruleElkNodeEdges[$current] )
                    // InternalElkGraphJson.g:249:4: ruleKeyEdges otherlv_12= ':' this_ElkNodeEdges_13= ruleElkNodeEdges[$current]
                    {

                    				newCompositeNode(grammarAccess.getNodeElementAccess().getKeyEdgesParserRuleCall_5_0());
                    			
                    pushFollow(FOLLOW_7);
                    ruleKeyEdges();

                    state._fsp--;


                    				afterParserOrEnumRuleCall();
                    			
                    otherlv_12=(Token)match(input,16,FOLLOW_8); 

                    				newLeafNode(otherlv_12, grammarAccess.getNodeElementAccess().getColonKeyword_5_1());
                    			

                    				if (current==null) {
                    					current = createModelElement(grammarAccess.getNodeElementRule());
                    				}
                    				newCompositeNode(grammarAccess.getNodeElementAccess().getElkNodeEdgesParserRuleCall_5_2());
                    			
                    pushFollow(FOLLOW_2);
                    this_ElkNodeEdges_13=ruleElkNodeEdges(current);

                    state._fsp--;


                    				current = this_ElkNodeEdges_13;
                    				afterParserOrEnumRuleCall();
                    			

                    }


                    }
                    break;
                case 7 :
                    // InternalElkGraphJson.g:273:3: ( ruleKeyLayoutOptions otherlv_15= ':' this_ElkGraphElementProperties_16= ruleElkGraphElementProperties[$current] )
                    {
                    // InternalElkGraphJson.g:273:3: ( ruleKeyLayoutOptions otherlv_15= ':' this_ElkGraphElementProperties_16= ruleElkGraphElementProperties[$current] )
                    // InternalElkGraphJson.g:274:4: ruleKeyLayoutOptions otherlv_15= ':' this_ElkGraphElementProperties_16= ruleElkGraphElementProperties[$current]
                    {

                    				newCompositeNode(grammarAccess.getNodeElementAccess().getKeyLayoutOptionsParserRuleCall_6_0());
                    			
                    pushFollow(FOLLOW_7);
                    ruleKeyLayoutOptions();

                    state._fsp--;


                    				afterParserOrEnumRuleCall();
                    			
                    otherlv_15=(Token)match(input,16,FOLLOW_9); 

                    				newLeafNode(otherlv_15, grammarAccess.getNodeElementAccess().getColonKeyword_6_1());
                    			

                    				if (current==null) {
                    					current = createModelElement(grammarAccess.getNodeElementRule());
                    				}
                    				newCompositeNode(grammarAccess.getNodeElementAccess().getElkGraphElementPropertiesParserRuleCall_6_2());
                    			
                    pushFollow(FOLLOW_2);
                    this_ElkGraphElementProperties_16=ruleElkGraphElementProperties(current);

                    state._fsp--;


                    				current = this_ElkGraphElementProperties_16;
                    				afterParserOrEnumRuleCall();
                    			

                    }


                    }
                    break;
                case 8 :
                    // InternalElkGraphJson.g:298:3: ruleJsonMember
                    {

                    			newCompositeNode(grammarAccess.getNodeElementAccess().getJsonMemberParserRuleCall_7());
                    		
                    pushFollow(FOLLOW_2);
                    ruleJsonMember();

                    state._fsp--;


                    			afterParserOrEnumRuleCall();
                    		

                    }
                    break;

            }


            }


            	leaveRule();

        }

            catch (RecognitionException re) {
                recover(input,re);
                appendSkippedTokens();
            }
        finally {
        }
        return current;
    }
    // $ANTLR end "ruleNodeElement"


    // $ANTLR start "entryRuleElkPort"
    // InternalElkGraphJson.g:309:1: entryRuleElkPort returns [EObject current=null] : iv_ruleElkPort= ruleElkPort EOF ;
    public final EObject entryRuleElkPort() throws RecognitionException {
        EObject current = null;

        EObject iv_ruleElkPort = null;


        try {
            // InternalElkGraphJson.g:309:48: (iv_ruleElkPort= ruleElkPort EOF )
            // InternalElkGraphJson.g:310:2: iv_ruleElkPort= ruleElkPort EOF
            {
             newCompositeNode(grammarAccess.getElkPortRule()); 
            pushFollow(FOLLOW_1);
            iv_ruleElkPort=ruleElkPort();

            state._fsp--;

             current =iv_ruleElkPort; 
            match(input,EOF,FOLLOW_2); 

            }

        }

            catch (RecognitionException re) {
                recover(input,re);
                appendSkippedTokens();
            }
        finally {
        }
        return current;
    }
    // $ANTLR end "entryRuleElkPort"


    // $ANTLR start "ruleElkPort"
    // InternalElkGraphJson.g:316:1: ruleElkPort returns [EObject current=null] : (otherlv_0= '{' this_PortElement_1= rulePortElement[$current] (otherlv_2= ',' this_PortElement_3= rulePortElement[$current] )* (otherlv_4= ',' )? otherlv_5= '}' ) ;
    public final EObject ruleElkPort() throws RecognitionException {
        EObject current = null;

        Token otherlv_0=null;
        Token otherlv_2=null;
        Token otherlv_4=null;
        Token otherlv_5=null;
        EObject this_PortElement_1 = null;

        EObject this_PortElement_3 = null;



        	enterRule();

        try {
            // InternalElkGraphJson.g:322:2: ( (otherlv_0= '{' this_PortElement_1= rulePortElement[$current] (otherlv_2= ',' this_PortElement_3= rulePortElement[$current] )* (otherlv_4= ',' )? otherlv_5= '}' ) )
            // InternalElkGraphJson.g:323:2: (otherlv_0= '{' this_PortElement_1= rulePortElement[$current] (otherlv_2= ',' this_PortElement_3= rulePortElement[$current] )* (otherlv_4= ',' )? otherlv_5= '}' )
            {
            // InternalElkGraphJson.g:323:2: (otherlv_0= '{' this_PortElement_1= rulePortElement[$current] (otherlv_2= ',' this_PortElement_3= rulePortElement[$current] )* (otherlv_4= ',' )? otherlv_5= '}' )
            // InternalElkGraphJson.g:324:3: otherlv_0= '{' this_PortElement_1= rulePortElement[$current] (otherlv_2= ',' this_PortElement_3= rulePortElement[$current] )* (otherlv_4= ',' )? otherlv_5= '}'
            {
            otherlv_0=(Token)match(input,13,FOLLOW_5); 

            			newLeafNode(otherlv_0, grammarAccess.getElkPortAccess().getLeftCurlyBracketKeyword_0());
            		

            			if (current==null) {
            				current = createModelElement(grammarAccess.getElkPortRule());
            			}
            			newCompositeNode(grammarAccess.getElkPortAccess().getPortElementParserRuleCall_1());
            		
            pushFollow(FOLLOW_4);
            this_PortElement_1=rulePortElement(current);

            state._fsp--;


            			current = this_PortElement_1;
            			afterParserOrEnumRuleCall();
            		
            // InternalElkGraphJson.g:339:3: (otherlv_2= ',' this_PortElement_3= rulePortElement[$current] )*
            loop5:
            do {
                int alt5=2;
                int LA5_0 = input.LA(1);

                if ( (LA5_0==14) ) {
                    int LA5_1 = input.LA(2);

                    if ( ((LA5_1>=RULE_STRING && LA5_1<=RULE_ID)||(LA5_1>=28 && LA5_1<=30)||(LA5_1>=34 && LA5_1<=54)) ) {
                        alt5=1;
                    }


                }


                switch (alt5) {
            	case 1 :
            	    // InternalElkGraphJson.g:340:4: otherlv_2= ',' this_PortElement_3= rulePortElement[$current]
            	    {
            	    otherlv_2=(Token)match(input,14,FOLLOW_5); 

            	    				newLeafNode(otherlv_2, grammarAccess.getElkPortAccess().getCommaKeyword_2_0());
            	    			

            	    				if (current==null) {
            	    					current = createModelElement(grammarAccess.getElkPortRule());
            	    				}
            	    				newCompositeNode(grammarAccess.getElkPortAccess().getPortElementParserRuleCall_2_1());
            	    			
            	    pushFollow(FOLLOW_4);
            	    this_PortElement_3=rulePortElement(current);

            	    state._fsp--;


            	    				current = this_PortElement_3;
            	    				afterParserOrEnumRuleCall();
            	    			

            	    }
            	    break;

            	default :
            	    break loop5;
                }
            } while (true);

            // InternalElkGraphJson.g:356:3: (otherlv_4= ',' )?
            int alt6=2;
            int LA6_0 = input.LA(1);

            if ( (LA6_0==14) ) {
                alt6=1;
            }
            switch (alt6) {
                case 1 :
                    // InternalElkGraphJson.g:357:4: otherlv_4= ','
                    {
                    otherlv_4=(Token)match(input,14,FOLLOW_6); 

                    				newLeafNode(otherlv_4, grammarAccess.getElkPortAccess().getCommaKeyword_3());
                    			

                    }
                    break;

            }

            otherlv_5=(Token)match(input,15,FOLLOW_2); 

            			newLeafNode(otherlv_5, grammarAccess.getElkPortAccess().getRightCurlyBracketKeyword_4());
            		

            }


            }


            	leaveRule();

        }

            catch (RecognitionException re) {
                recover(input,re);
                appendSkippedTokens();
            }
        finally {
        }
        return current;
    }
    // $ANTLR end "ruleElkPort"


    // $ANTLR start "rulePortElement"
    // InternalElkGraphJson.g:371:1: rulePortElement[EObject in_current] returns [EObject current=in_current] : (this_ElkId_0= ruleElkId[$current] | this_ShapeElement_1= ruleShapeElement[$current] | ( ruleKeyLabels otherlv_3= ':' this_ElkGraphElementLabels_4= ruleElkGraphElementLabels[$current] ) | ( ruleKeyLayoutOptions otherlv_6= ':' this_ElkGraphElementProperties_7= ruleElkGraphElementProperties[$current] ) | ruleJsonMember ) ;
    public final EObject rulePortElement(EObject in_current) throws RecognitionException {
        EObject current = in_current;

        Token otherlv_3=null;
        Token otherlv_6=null;
        EObject this_ElkId_0 = null;

        EObject this_ShapeElement_1 = null;

        EObject this_ElkGraphElementLabels_4 = null;

        EObject this_ElkGraphElementProperties_7 = null;



        	enterRule();

        try {
            // InternalElkGraphJson.g:377:2: ( (this_ElkId_0= ruleElkId[$current] | this_ShapeElement_1= ruleShapeElement[$current] | ( ruleKeyLabels otherlv_3= ':' this_ElkGraphElementLabels_4= ruleElkGraphElementLabels[$current] ) | ( ruleKeyLayoutOptions otherlv_6= ':' this_ElkGraphElementProperties_7= ruleElkGraphElementProperties[$current] ) | ruleJsonMember ) )
            // InternalElkGraphJson.g:378:2: (this_ElkId_0= ruleElkId[$current] | this_ShapeElement_1= ruleShapeElement[$current] | ( ruleKeyLabels otherlv_3= ':' this_ElkGraphElementLabels_4= ruleElkGraphElementLabels[$current] ) | ( ruleKeyLayoutOptions otherlv_6= ':' this_ElkGraphElementProperties_7= ruleElkGraphElementProperties[$current] ) | ruleJsonMember )
            {
            // InternalElkGraphJson.g:378:2: (this_ElkId_0= ruleElkId[$current] | this_ShapeElement_1= ruleShapeElement[$current] | ( ruleKeyLabels otherlv_3= ':' this_ElkGraphElementLabels_4= ruleElkGraphElementLabels[$current] ) | ( ruleKeyLayoutOptions otherlv_6= ':' this_ElkGraphElementProperties_7= ruleElkGraphElementProperties[$current] ) | ruleJsonMember )
            int alt7=5;
            switch ( input.LA(1) ) {
            case 40:
            case 41:
            case 42:
                {
                alt7=1;
                }
                break;
            case 43:
            case 44:
            case 45:
            case 46:
            case 47:
            case 48:
            case 49:
            case 50:
            case 51:
            case 52:
            case 53:
            case 54:
                {
                alt7=2;
                }
                break;
            case 28:
            case 29:
            case 30:
                {
                alt7=3;
                }
                break;
            case 34:
            case 35:
            case 36:
            case 37:
            case 38:
            case 39:
                {
                alt7=4;
                }
                break;
            case RULE_STRING:
            case RULE_ID:
                {
                alt7=5;
                }
                break;
            default:
                NoViableAltException nvae =
                    new NoViableAltException("", 7, 0, input);

                throw nvae;
            }

            switch (alt7) {
                case 1 :
                    // InternalElkGraphJson.g:379:3: this_ElkId_0= ruleElkId[$current]
                    {

                    			if (current==null) {
                    				current = createModelElement(grammarAccess.getPortElementRule());
                    			}
                    			newCompositeNode(grammarAccess.getPortElementAccess().getElkIdParserRuleCall_0());
                    		
                    pushFollow(FOLLOW_2);
                    this_ElkId_0=ruleElkId(current);

                    state._fsp--;


                    			current = this_ElkId_0;
                    			afterParserOrEnumRuleCall();
                    		

                    }
                    break;
                case 2 :
                    // InternalElkGraphJson.g:391:3: this_ShapeElement_1= ruleShapeElement[$current]
                    {

                    			if (current==null) {
                    				current = createModelElement(grammarAccess.getPortElementRule());
                    			}
                    			newCompositeNode(grammarAccess.getPortElementAccess().getShapeElementParserRuleCall_1());
                    		
                    pushFollow(FOLLOW_2);
                    this_ShapeElement_1=ruleShapeElement(current);

                    state._fsp--;


                    			current = this_ShapeElement_1;
                    			afterParserOrEnumRuleCall();
                    		

                    }
                    break;
                case 3 :
                    // InternalElkGraphJson.g:403:3: ( ruleKeyLabels otherlv_3= ':' this_ElkGraphElementLabels_4= ruleElkGraphElementLabels[$current] )
                    {
                    // InternalElkGraphJson.g:403:3: ( ruleKeyLabels otherlv_3= ':' this_ElkGraphElementLabels_4= ruleElkGraphElementLabels[$current] )
                    // InternalElkGraphJson.g:404:4: ruleKeyLabels otherlv_3= ':' this_ElkGraphElementLabels_4= ruleElkGraphElementLabels[$current]
                    {

                    				newCompositeNode(grammarAccess.getPortElementAccess().getKeyLabelsParserRuleCall_2_0());
                    			
                    pushFollow(FOLLOW_7);
                    ruleKeyLabels();

                    state._fsp--;


                    				afterParserOrEnumRuleCall();
                    			
                    otherlv_3=(Token)match(input,16,FOLLOW_8); 

                    				newLeafNode(otherlv_3, grammarAccess.getPortElementAccess().getColonKeyword_2_1());
                    			

                    				if (current==null) {
                    					current = createModelElement(grammarAccess.getPortElementRule());
                    				}
                    				newCompositeNode(grammarAccess.getPortElementAccess().getElkGraphElementLabelsParserRuleCall_2_2());
                    			
                    pushFollow(FOLLOW_2);
                    this_ElkGraphElementLabels_4=ruleElkGraphElementLabels(current);

                    state._fsp--;


                    				current = this_ElkGraphElementLabels_4;
                    				afterParserOrEnumRuleCall();
                    			

                    }


                    }
                    break;
                case 4 :
                    // InternalElkGraphJson.g:428:3: ( ruleKeyLayoutOptions otherlv_6= ':' this_ElkGraphElementProperties_7= ruleElkGraphElementProperties[$current] )
                    {
                    // InternalElkGraphJson.g:428:3: ( ruleKeyLayoutOptions otherlv_6= ':' this_ElkGraphElementProperties_7= ruleElkGraphElementProperties[$current] )
                    // InternalElkGraphJson.g:429:4: ruleKeyLayoutOptions otherlv_6= ':' this_ElkGraphElementProperties_7= ruleElkGraphElementProperties[$current]
                    {

                    				newCompositeNode(grammarAccess.getPortElementAccess().getKeyLayoutOptionsParserRuleCall_3_0());
                    			
                    pushFollow(FOLLOW_7);
                    ruleKeyLayoutOptions();

                    state._fsp--;


                    				afterParserOrEnumRuleCall();
                    			
                    otherlv_6=(Token)match(input,16,FOLLOW_9); 

                    				newLeafNode(otherlv_6, grammarAccess.getPortElementAccess().getColonKeyword_3_1());
                    			

                    				if (current==null) {
                    					current = createModelElement(grammarAccess.getPortElementRule());
                    				}
                    				newCompositeNode(grammarAccess.getPortElementAccess().getElkGraphElementPropertiesParserRuleCall_3_2());
                    			
                    pushFollow(FOLLOW_2);
                    this_ElkGraphElementProperties_7=ruleElkGraphElementProperties(current);

                    state._fsp--;


                    				current = this_ElkGraphElementProperties_7;
                    				afterParserOrEnumRuleCall();
                    			

                    }


                    }
                    break;
                case 5 :
                    // InternalElkGraphJson.g:453:3: ruleJsonMember
                    {

                    			newCompositeNode(grammarAccess.getPortElementAccess().getJsonMemberParserRuleCall_4());
                    		
                    pushFollow(FOLLOW_2);
                    ruleJsonMember();

                    state._fsp--;


                    			afterParserOrEnumRuleCall();
                    		

                    }
                    break;

            }


            }


            	leaveRule();

        }

            catch (RecognitionException re) {
                recover(input,re);
                appendSkippedTokens();
            }
        finally {
        }
        return current;
    }
    // $ANTLR end "rulePortElement"


    // $ANTLR start "entryRuleElkLabel"
    // InternalElkGraphJson.g:464:1: entryRuleElkLabel returns [EObject current=null] : iv_ruleElkLabel= ruleElkLabel EOF ;
    public final EObject entryRuleElkLabel() throws RecognitionException {
        EObject current = null;

        EObject iv_ruleElkLabel = null;


        try {
            // InternalElkGraphJson.g:464:49: (iv_ruleElkLabel= ruleElkLabel EOF )
            // InternalElkGraphJson.g:465:2: iv_ruleElkLabel= ruleElkLabel EOF
            {
             newCompositeNode(grammarAccess.getElkLabelRule()); 
            pushFollow(FOLLOW_1);
            iv_ruleElkLabel=ruleElkLabel();

            state._fsp--;

             current =iv_ruleElkLabel; 
            match(input,EOF,FOLLOW_2); 

            }

        }

            catch (RecognitionException re) {
                recover(input,re);
                appendSkippedTokens();
            }
        finally {
        }
        return current;
    }
    // $ANTLR end "entryRuleElkLabel"


    // $ANTLR start "ruleElkLabel"
    // InternalElkGraphJson.g:471:1: ruleElkLabel returns [EObject current=null] : (otherlv_0= '{' this_LabelElement_1= ruleLabelElement[$current] (otherlv_2= ',' this_LabelElement_3= ruleLabelElement[$current] )* (otherlv_4= ',' )? otherlv_5= '}' ) ;
    public final EObject ruleElkLabel() throws RecognitionException {
        EObject current = null;

        Token otherlv_0=null;
        Token otherlv_2=null;
        Token otherlv_4=null;
        Token otherlv_5=null;
        EObject this_LabelElement_1 = null;

        EObject this_LabelElement_3 = null;



        	enterRule();

        try {
            // InternalElkGraphJson.g:477:2: ( (otherlv_0= '{' this_LabelElement_1= ruleLabelElement[$current] (otherlv_2= ',' this_LabelElement_3= ruleLabelElement[$current] )* (otherlv_4= ',' )? otherlv_5= '}' ) )
            // InternalElkGraphJson.g:478:2: (otherlv_0= '{' this_LabelElement_1= ruleLabelElement[$current] (otherlv_2= ',' this_LabelElement_3= ruleLabelElement[$current] )* (otherlv_4= ',' )? otherlv_5= '}' )
            {
            // InternalElkGraphJson.g:478:2: (otherlv_0= '{' this_LabelElement_1= ruleLabelElement[$current] (otherlv_2= ',' this_LabelElement_3= ruleLabelElement[$current] )* (otherlv_4= ',' )? otherlv_5= '}' )
            // InternalElkGraphJson.g:479:3: otherlv_0= '{' this_LabelElement_1= ruleLabelElement[$current] (otherlv_2= ',' this_LabelElement_3= ruleLabelElement[$current] )* (otherlv_4= ',' )? otherlv_5= '}'
            {
            otherlv_0=(Token)match(input,13,FOLLOW_10); 

            			newLeafNode(otherlv_0, grammarAccess.getElkLabelAccess().getLeftCurlyBracketKeyword_0());
            		

            			if (current==null) {
            				current = createModelElement(grammarAccess.getElkLabelRule());
            			}
            			newCompositeNode(grammarAccess.getElkLabelAccess().getLabelElementParserRuleCall_1());
            		
            pushFollow(FOLLOW_4);
            this_LabelElement_1=ruleLabelElement(current);

            state._fsp--;


            			current = this_LabelElement_1;
            			afterParserOrEnumRuleCall();
            		
            // InternalElkGraphJson.g:494:3: (otherlv_2= ',' this_LabelElement_3= ruleLabelElement[$current] )*
            loop8:
            do {
                int alt8=2;
                int LA8_0 = input.LA(1);

                if ( (LA8_0==14) ) {
                    int LA8_1 = input.LA(2);

                    if ( ((LA8_1>=RULE_STRING && LA8_1<=RULE_ID)||(LA8_1>=28 && LA8_1<=30)||(LA8_1>=34 && LA8_1<=54)||(LA8_1>=61 && LA8_1<=63)) ) {
                        alt8=1;
                    }


                }


                switch (alt8) {
            	case 1 :
            	    // InternalElkGraphJson.g:495:4: otherlv_2= ',' this_LabelElement_3= ruleLabelElement[$current]
            	    {
            	    otherlv_2=(Token)match(input,14,FOLLOW_10); 

            	    				newLeafNode(otherlv_2, grammarAccess.getElkLabelAccess().getCommaKeyword_2_0());
            	    			

            	    				if (current==null) {
            	    					current = createModelElement(grammarAccess.getElkLabelRule());
            	    				}
            	    				newCompositeNode(grammarAccess.getElkLabelAccess().getLabelElementParserRuleCall_2_1());
            	    			
            	    pushFollow(FOLLOW_4);
            	    this_LabelElement_3=ruleLabelElement(current);

            	    state._fsp--;


            	    				current = this_LabelElement_3;
            	    				afterParserOrEnumRuleCall();
            	    			

            	    }
            	    break;

            	default :
            	    break loop8;
                }
            } while (true);

            // InternalElkGraphJson.g:511:3: (otherlv_4= ',' )?
            int alt9=2;
            int LA9_0 = input.LA(1);

            if ( (LA9_0==14) ) {
                alt9=1;
            }
            switch (alt9) {
                case 1 :
                    // InternalElkGraphJson.g:512:4: otherlv_4= ','
                    {
                    otherlv_4=(Token)match(input,14,FOLLOW_6); 

                    				newLeafNode(otherlv_4, grammarAccess.getElkLabelAccess().getCommaKeyword_3());
                    			

                    }
                    break;

            }

            otherlv_5=(Token)match(input,15,FOLLOW_2); 

            			newLeafNode(otherlv_5, grammarAccess.getElkLabelAccess().getRightCurlyBracketKeyword_4());
            		

            }


            }


            	leaveRule();

        }

            catch (RecognitionException re) {
                recover(input,re);
                appendSkippedTokens();
            }
        finally {
        }
        return current;
    }
    // $ANTLR end "ruleElkLabel"


    // $ANTLR start "ruleLabelElement"
    // InternalElkGraphJson.g:526:1: ruleLabelElement[EObject in_current] returns [EObject current=in_current] : (this_ElkId_0= ruleElkId[$current] | this_ShapeElement_1= ruleShapeElement[$current] | ( ruleKeyText otherlv_3= ':' ( (lv_text_4_0= RULE_STRING ) ) ) | ( ruleKeyLabels otherlv_6= ':' this_ElkGraphElementLabels_7= ruleElkGraphElementLabels[$current] ) | ( ruleKeyLayoutOptions otherlv_9= ':' this_ElkGraphElementProperties_10= ruleElkGraphElementProperties[$current] ) | ruleJsonMember ) ;
    public final EObject ruleLabelElement(EObject in_current) throws RecognitionException {
        EObject current = in_current;

        Token otherlv_3=null;
        Token lv_text_4_0=null;
        Token otherlv_6=null;
        Token otherlv_9=null;
        EObject this_ElkId_0 = null;

        EObject this_ShapeElement_1 = null;

        EObject this_ElkGraphElementLabels_7 = null;

        EObject this_ElkGraphElementProperties_10 = null;



        	enterRule();

        try {
            // InternalElkGraphJson.g:532:2: ( (this_ElkId_0= ruleElkId[$current] | this_ShapeElement_1= ruleShapeElement[$current] | ( ruleKeyText otherlv_3= ':' ( (lv_text_4_0= RULE_STRING ) ) ) | ( ruleKeyLabels otherlv_6= ':' this_ElkGraphElementLabels_7= ruleElkGraphElementLabels[$current] ) | ( ruleKeyLayoutOptions otherlv_9= ':' this_ElkGraphElementProperties_10= ruleElkGraphElementProperties[$current] ) | ruleJsonMember ) )
            // InternalElkGraphJson.g:533:2: (this_ElkId_0= ruleElkId[$current] | this_ShapeElement_1= ruleShapeElement[$current] | ( ruleKeyText otherlv_3= ':' ( (lv_text_4_0= RULE_STRING ) ) ) | ( ruleKeyLabels otherlv_6= ':' this_ElkGraphElementLabels_7= ruleElkGraphElementLabels[$current] ) | ( ruleKeyLayoutOptions otherlv_9= ':' this_ElkGraphElementProperties_10= ruleElkGraphElementProperties[$current] ) | ruleJsonMember )
            {
            // InternalElkGraphJson.g:533:2: (this_ElkId_0= ruleElkId[$current] | this_ShapeElement_1= ruleShapeElement[$current] | ( ruleKeyText otherlv_3= ':' ( (lv_text_4_0= RULE_STRING ) ) ) | ( ruleKeyLabels otherlv_6= ':' this_ElkGraphElementLabels_7= ruleElkGraphElementLabels[$current] ) | ( ruleKeyLayoutOptions otherlv_9= ':' this_ElkGraphElementProperties_10= ruleElkGraphElementProperties[$current] ) | ruleJsonMember )
            int alt10=6;
            switch ( input.LA(1) ) {
            case 40:
            case 41:
            case 42:
                {
                alt10=1;
                }
                break;
            case 43:
            case 44:
            case 45:
            case 46:
            case 47:
            case 48:
            case 49:
            case 50:
            case 51:
            case 52:
            case 53:
            case 54:
                {
                alt10=2;
                }
                break;
            case 61:
            case 62:
            case 63:
                {
                alt10=3;
                }
                break;
            case 28:
            case 29:
            case 30:
                {
                alt10=4;
                }
                break;
            case 34:
            case 35:
            case 36:
            case 37:
            case 38:
            case 39:
                {
                alt10=5;
                }
                break;
            case RULE_STRING:
            case RULE_ID:
                {
                alt10=6;
                }
                break;
            default:
                NoViableAltException nvae =
                    new NoViableAltException("", 10, 0, input);

                throw nvae;
            }

            switch (alt10) {
                case 1 :
                    // InternalElkGraphJson.g:534:3: this_ElkId_0= ruleElkId[$current]
                    {

                    			if (current==null) {
                    				current = createModelElement(grammarAccess.getLabelElementRule());
                    			}
                    			newCompositeNode(grammarAccess.getLabelElementAccess().getElkIdParserRuleCall_0());
                    		
                    pushFollow(FOLLOW_2);
                    this_ElkId_0=ruleElkId(current);

                    state._fsp--;


                    			current = this_ElkId_0;
                    			afterParserOrEnumRuleCall();
                    		

                    }
                    break;
                case 2 :
                    // InternalElkGraphJson.g:546:3: this_ShapeElement_1= ruleShapeElement[$current]
                    {

                    			if (current==null) {
                    				current = createModelElement(grammarAccess.getLabelElementRule());
                    			}
                    			newCompositeNode(grammarAccess.getLabelElementAccess().getShapeElementParserRuleCall_1());
                    		
                    pushFollow(FOLLOW_2);
                    this_ShapeElement_1=ruleShapeElement(current);

                    state._fsp--;


                    			current = this_ShapeElement_1;
                    			afterParserOrEnumRuleCall();
                    		

                    }
                    break;
                case 3 :
                    // InternalElkGraphJson.g:558:3: ( ruleKeyText otherlv_3= ':' ( (lv_text_4_0= RULE_STRING ) ) )
                    {
                    // InternalElkGraphJson.g:558:3: ( ruleKeyText otherlv_3= ':' ( (lv_text_4_0= RULE_STRING ) ) )
                    // InternalElkGraphJson.g:559:4: ruleKeyText otherlv_3= ':' ( (lv_text_4_0= RULE_STRING ) )
                    {

                    				newCompositeNode(grammarAccess.getLabelElementAccess().getKeyTextParserRuleCall_2_0());
                    			
                    pushFollow(FOLLOW_7);
                    ruleKeyText();

                    state._fsp--;


                    				afterParserOrEnumRuleCall();
                    			
                    otherlv_3=(Token)match(input,16,FOLLOW_11); 

                    				newLeafNode(otherlv_3, grammarAccess.getLabelElementAccess().getColonKeyword_2_1());
                    			
                    // InternalElkGraphJson.g:570:4: ( (lv_text_4_0= RULE_STRING ) )
                    // InternalElkGraphJson.g:571:5: (lv_text_4_0= RULE_STRING )
                    {
                    // InternalElkGraphJson.g:571:5: (lv_text_4_0= RULE_STRING )
                    // InternalElkGraphJson.g:572:6: lv_text_4_0= RULE_STRING
                    {
                    lv_text_4_0=(Token)match(input,RULE_STRING,FOLLOW_2); 

                    						newLeafNode(lv_text_4_0, grammarAccess.getLabelElementAccess().getTextSTRINGTerminalRuleCall_2_2_0());
                    					

                    						if (current==null) {
                    							current = createModelElement(grammarAccess.getLabelElementRule());
                    						}
                    						setWithLastConsumed(
                    							current,
                    							"text",
                    							lv_text_4_0,
                    							"org.eclipse.xtext.common.Terminals.STRING");
                    					

                    }


                    }


                    }


                    }
                    break;
                case 4 :
                    // InternalElkGraphJson.g:590:3: ( ruleKeyLabels otherlv_6= ':' this_ElkGraphElementLabels_7= ruleElkGraphElementLabels[$current] )
                    {
                    // InternalElkGraphJson.g:590:3: ( ruleKeyLabels otherlv_6= ':' this_ElkGraphElementLabels_7= ruleElkGraphElementLabels[$current] )
                    // InternalElkGraphJson.g:591:4: ruleKeyLabels otherlv_6= ':' this_ElkGraphElementLabels_7= ruleElkGraphElementLabels[$current]
                    {

                    				newCompositeNode(grammarAccess.getLabelElementAccess().getKeyLabelsParserRuleCall_3_0());
                    			
                    pushFollow(FOLLOW_7);
                    ruleKeyLabels();

                    state._fsp--;


                    				afterParserOrEnumRuleCall();
                    			
                    otherlv_6=(Token)match(input,16,FOLLOW_8); 

                    				newLeafNode(otherlv_6, grammarAccess.getLabelElementAccess().getColonKeyword_3_1());
                    			

                    				if (current==null) {
                    					current = createModelElement(grammarAccess.getLabelElementRule());
                    				}
                    				newCompositeNode(grammarAccess.getLabelElementAccess().getElkGraphElementLabelsParserRuleCall_3_2());
                    			
                    pushFollow(FOLLOW_2);
                    this_ElkGraphElementLabels_7=ruleElkGraphElementLabels(current);

                    state._fsp--;


                    				current = this_ElkGraphElementLabels_7;
                    				afterParserOrEnumRuleCall();
                    			

                    }


                    }
                    break;
                case 5 :
                    // InternalElkGraphJson.g:615:3: ( ruleKeyLayoutOptions otherlv_9= ':' this_ElkGraphElementProperties_10= ruleElkGraphElementProperties[$current] )
                    {
                    // InternalElkGraphJson.g:615:3: ( ruleKeyLayoutOptions otherlv_9= ':' this_ElkGraphElementProperties_10= ruleElkGraphElementProperties[$current] )
                    // InternalElkGraphJson.g:616:4: ruleKeyLayoutOptions otherlv_9= ':' this_ElkGraphElementProperties_10= ruleElkGraphElementProperties[$current]
                    {

                    				newCompositeNode(grammarAccess.getLabelElementAccess().getKeyLayoutOptionsParserRuleCall_4_0());
                    			
                    pushFollow(FOLLOW_7);
                    ruleKeyLayoutOptions();

                    state._fsp--;


                    				afterParserOrEnumRuleCall();
                    			
                    otherlv_9=(Token)match(input,16,FOLLOW_9); 

                    				newLeafNode(otherlv_9, grammarAccess.getLabelElementAccess().getColonKeyword_4_1());
                    			

                    				if (current==null) {
                    					current = createModelElement(grammarAccess.getLabelElementRule());
                    				}
                    				newCompositeNode(grammarAccess.getLabelElementAccess().getElkGraphElementPropertiesParserRuleCall_4_2());
                    			
                    pushFollow(FOLLOW_2);
                    this_ElkGraphElementProperties_10=ruleElkGraphElementProperties(current);

                    state._fsp--;


                    				current = this_ElkGraphElementProperties_10;
                    				afterParserOrEnumRuleCall();
                    			

                    }


                    }
                    break;
                case 6 :
                    // InternalElkGraphJson.g:640:3: ruleJsonMember
                    {

                    			newCompositeNode(grammarAccess.getLabelElementAccess().getJsonMemberParserRuleCall_5());
                    		
                    pushFollow(FOLLOW_2);
                    ruleJsonMember();

                    state._fsp--;


                    			afterParserOrEnumRuleCall();
                    		

                    }
                    break;

            }


            }


            	leaveRule();

        }

            catch (RecognitionException re) {
                recover(input,re);
                appendSkippedTokens();
            }
        finally {
        }
        return current;
    }
    // $ANTLR end "ruleLabelElement"


    // $ANTLR start "entryRuleElkEdge"
    // InternalElkGraphJson.g:651:1: entryRuleElkEdge returns [EObject current=null] : iv_ruleElkEdge= ruleElkEdge EOF ;
    public final EObject entryRuleElkEdge() throws RecognitionException {
        EObject current = null;

        EObject iv_ruleElkEdge = null;


        try {
            // InternalElkGraphJson.g:651:48: (iv_ruleElkEdge= ruleElkEdge EOF )
            // InternalElkGraphJson.g:652:2: iv_ruleElkEdge= ruleElkEdge EOF
            {
             newCompositeNode(grammarAccess.getElkEdgeRule()); 
            pushFollow(FOLLOW_1);
            iv_ruleElkEdge=ruleElkEdge();

            state._fsp--;

             current =iv_ruleElkEdge; 
            match(input,EOF,FOLLOW_2); 

            }

        }

            catch (RecognitionException re) {
                recover(input,re);
                appendSkippedTokens();
            }
        finally {
        }
        return current;
    }
    // $ANTLR end "entryRuleElkEdge"


    // $ANTLR start "ruleElkEdge"
    // InternalElkGraphJson.g:658:1: ruleElkEdge returns [EObject current=null] : (otherlv_0= '{' this_EdgeElement_1= ruleEdgeElement[$current] (otherlv_2= ',' this_EdgeElement_3= ruleEdgeElement[$current] )* (otherlv_4= ',' )? otherlv_5= '}' ) ;
    public final EObject ruleElkEdge() throws RecognitionException {
        EObject current = null;

        Token otherlv_0=null;
        Token otherlv_2=null;
        Token otherlv_4=null;
        Token otherlv_5=null;
        EObject this_EdgeElement_1 = null;

        EObject this_EdgeElement_3 = null;



        	enterRule();

        try {
            // InternalElkGraphJson.g:664:2: ( (otherlv_0= '{' this_EdgeElement_1= ruleEdgeElement[$current] (otherlv_2= ',' this_EdgeElement_3= ruleEdgeElement[$current] )* (otherlv_4= ',' )? otherlv_5= '}' ) )
            // InternalElkGraphJson.g:665:2: (otherlv_0= '{' this_EdgeElement_1= ruleEdgeElement[$current] (otherlv_2= ',' this_EdgeElement_3= ruleEdgeElement[$current] )* (otherlv_4= ',' )? otherlv_5= '}' )
            {
            // InternalElkGraphJson.g:665:2: (otherlv_0= '{' this_EdgeElement_1= ruleEdgeElement[$current] (otherlv_2= ',' this_EdgeElement_3= ruleEdgeElement[$current] )* (otherlv_4= ',' )? otherlv_5= '}' )
            // InternalElkGraphJson.g:666:3: otherlv_0= '{' this_EdgeElement_1= ruleEdgeElement[$current] (otherlv_2= ',' this_EdgeElement_3= ruleEdgeElement[$current] )* (otherlv_4= ',' )? otherlv_5= '}'
            {
            otherlv_0=(Token)match(input,13,FOLLOW_12); 

            			newLeafNode(otherlv_0, grammarAccess.getElkEdgeAccess().getLeftCurlyBracketKeyword_0());
            		

            			if (current==null) {
            				current = createModelElement(grammarAccess.getElkEdgeRule());
            			}
            			newCompositeNode(grammarAccess.getElkEdgeAccess().getEdgeElementParserRuleCall_1());
            		
            pushFollow(FOLLOW_4);
            this_EdgeElement_1=ruleEdgeElement(current);

            state._fsp--;


            			current = this_EdgeElement_1;
            			afterParserOrEnumRuleCall();
            		
            // InternalElkGraphJson.g:681:3: (otherlv_2= ',' this_EdgeElement_3= ruleEdgeElement[$current] )*
            loop11:
            do {
                int alt11=2;
                int LA11_0 = input.LA(1);

                if ( (LA11_0==14) ) {
                    int LA11_1 = input.LA(2);

                    if ( ((LA11_1>=RULE_STRING && LA11_1<=RULE_ID)||(LA11_1>=28 && LA11_1<=30)||(LA11_1>=34 && LA11_1<=42)||(LA11_1>=55 && LA11_1<=60)) ) {
                        alt11=1;
                    }


                }


                switch (alt11) {
            	case 1 :
            	    // InternalElkGraphJson.g:682:4: otherlv_2= ',' this_EdgeElement_3= ruleEdgeElement[$current]
            	    {
            	    otherlv_2=(Token)match(input,14,FOLLOW_12); 

            	    				newLeafNode(otherlv_2, grammarAccess.getElkEdgeAccess().getCommaKeyword_2_0());
            	    			

            	    				if (current==null) {
            	    					current = createModelElement(grammarAccess.getElkEdgeRule());
            	    				}
            	    				newCompositeNode(grammarAccess.getElkEdgeAccess().getEdgeElementParserRuleCall_2_1());
            	    			
            	    pushFollow(FOLLOW_4);
            	    this_EdgeElement_3=ruleEdgeElement(current);

            	    state._fsp--;


            	    				current = this_EdgeElement_3;
            	    				afterParserOrEnumRuleCall();
            	    			

            	    }
            	    break;

            	default :
            	    break loop11;
                }
            } while (true);

            // InternalElkGraphJson.g:698:3: (otherlv_4= ',' )?
            int alt12=2;
            int LA12_0 = input.LA(1);

            if ( (LA12_0==14) ) {
                alt12=1;
            }
            switch (alt12) {
                case 1 :
                    // InternalElkGraphJson.g:699:4: otherlv_4= ','
                    {
                    otherlv_4=(Token)match(input,14,FOLLOW_6); 

                    				newLeafNode(otherlv_4, grammarAccess.getElkEdgeAccess().getCommaKeyword_3());
                    			

                    }
                    break;

            }

            otherlv_5=(Token)match(input,15,FOLLOW_2); 

            			newLeafNode(otherlv_5, grammarAccess.getElkEdgeAccess().getRightCurlyBracketKeyword_4());
            		

            }


            }


            	leaveRule();

        }

            catch (RecognitionException re) {
                recover(input,re);
                appendSkippedTokens();
            }
        finally {
        }
        return current;
    }
    // $ANTLR end "ruleElkEdge"


    // $ANTLR start "ruleEdgeElement"
    // InternalElkGraphJson.g:713:1: ruleEdgeElement[EObject in_current] returns [EObject current=in_current] : (this_ElkId_0= ruleElkId[$current] | ( ruleKeySources otherlv_2= ':' this_ElkEdgeSources_3= ruleElkEdgeSources[$current] ) | ( ruleKeyTargets otherlv_5= ':' this_ElkEdgeTargets_6= ruleElkEdgeTargets[$current] ) | ( ruleKeyLabels otherlv_8= ':' this_ElkGraphElementLabels_9= ruleElkGraphElementLabels[$current] ) | ( ruleKeyLayoutOptions otherlv_11= ':' this_ElkGraphElementProperties_12= ruleElkGraphElementProperties[$current] ) | ruleJsonMember ) ;
    public final EObject ruleEdgeElement(EObject in_current) throws RecognitionException {
        EObject current = in_current;

        Token otherlv_2=null;
        Token otherlv_5=null;
        Token otherlv_8=null;
        Token otherlv_11=null;
        EObject this_ElkId_0 = null;

        EObject this_ElkEdgeSources_3 = null;

        EObject this_ElkEdgeTargets_6 = null;

        EObject this_ElkGraphElementLabels_9 = null;

        EObject this_ElkGraphElementProperties_12 = null;



        	enterRule();

        try {
            // InternalElkGraphJson.g:719:2: ( (this_ElkId_0= ruleElkId[$current] | ( ruleKeySources otherlv_2= ':' this_ElkEdgeSources_3= ruleElkEdgeSources[$current] ) | ( ruleKeyTargets otherlv_5= ':' this_ElkEdgeTargets_6= ruleElkEdgeTargets[$current] ) | ( ruleKeyLabels otherlv_8= ':' this_ElkGraphElementLabels_9= ruleElkGraphElementLabels[$current] ) | ( ruleKeyLayoutOptions otherlv_11= ':' this_ElkGraphElementProperties_12= ruleElkGraphElementProperties[$current] ) | ruleJsonMember ) )
            // InternalElkGraphJson.g:720:2: (this_ElkId_0= ruleElkId[$current] | ( ruleKeySources otherlv_2= ':' this_ElkEdgeSources_3= ruleElkEdgeSources[$current] ) | ( ruleKeyTargets otherlv_5= ':' this_ElkEdgeTargets_6= ruleElkEdgeTargets[$current] ) | ( ruleKeyLabels otherlv_8= ':' this_ElkGraphElementLabels_9= ruleElkGraphElementLabels[$current] ) | ( ruleKeyLayoutOptions otherlv_11= ':' this_ElkGraphElementProperties_12= ruleElkGraphElementProperties[$current] ) | ruleJsonMember )
            {
            // InternalElkGraphJson.g:720:2: (this_ElkId_0= ruleElkId[$current] | ( ruleKeySources otherlv_2= ':' this_ElkEdgeSources_3= ruleElkEdgeSources[$current] ) | ( ruleKeyTargets otherlv_5= ':' this_ElkEdgeTargets_6= ruleElkEdgeTargets[$current] ) | ( ruleKeyLabels otherlv_8= ':' this_ElkGraphElementLabels_9= ruleElkGraphElementLabels[$current] ) | ( ruleKeyLayoutOptions otherlv_11= ':' this_ElkGraphElementProperties_12= ruleElkGraphElementProperties[$current] ) | ruleJsonMember )
            int alt13=6;
            switch ( input.LA(1) ) {
            case 40:
            case 41:
            case 42:
                {
                alt13=1;
                }
                break;
            case 55:
            case 56:
            case 57:
                {
                alt13=2;
                }
                break;
            case 58:
            case 59:
            case 60:
                {
                alt13=3;
                }
                break;
            case 28:
            case 29:
            case 30:
                {
                alt13=4;
                }
                break;
            case 34:
            case 35:
            case 36:
            case 37:
            case 38:
            case 39:
                {
                alt13=5;
                }
                break;
            case RULE_STRING:
            case RULE_ID:
                {
                alt13=6;
                }
                break;
            default:
                NoViableAltException nvae =
                    new NoViableAltException("", 13, 0, input);

                throw nvae;
            }

            switch (alt13) {
                case 1 :
                    // InternalElkGraphJson.g:721:3: this_ElkId_0= ruleElkId[$current]
                    {

                    			if (current==null) {
                    				current = createModelElement(grammarAccess.getEdgeElementRule());
                    			}
                    			newCompositeNode(grammarAccess.getEdgeElementAccess().getElkIdParserRuleCall_0());
                    		
                    pushFollow(FOLLOW_2);
                    this_ElkId_0=ruleElkId(current);

                    state._fsp--;


                    			current = this_ElkId_0;
                    			afterParserOrEnumRuleCall();
                    		

                    }
                    break;
                case 2 :
                    // InternalElkGraphJson.g:733:3: ( ruleKeySources otherlv_2= ':' this_ElkEdgeSources_3= ruleElkEdgeSources[$current] )
                    {
                    // InternalElkGraphJson.g:733:3: ( ruleKeySources otherlv_2= ':' this_ElkEdgeSources_3= ruleElkEdgeSources[$current] )
                    // InternalElkGraphJson.g:734:4: ruleKeySources otherlv_2= ':' this_ElkEdgeSources_3= ruleElkEdgeSources[$current]
                    {

                    				newCompositeNode(grammarAccess.getEdgeElementAccess().getKeySourcesParserRuleCall_1_0());
                    			
                    pushFollow(FOLLOW_7);
                    ruleKeySources();

                    state._fsp--;


                    				afterParserOrEnumRuleCall();
                    			
                    otherlv_2=(Token)match(input,16,FOLLOW_8); 

                    				newLeafNode(otherlv_2, grammarAccess.getEdgeElementAccess().getColonKeyword_1_1());
                    			

                    				if (current==null) {
                    					current = createModelElement(grammarAccess.getEdgeElementRule());
                    				}
                    				newCompositeNode(grammarAccess.getEdgeElementAccess().getElkEdgeSourcesParserRuleCall_1_2());
                    			
                    pushFollow(FOLLOW_2);
                    this_ElkEdgeSources_3=ruleElkEdgeSources(current);

                    state._fsp--;


                    				current = this_ElkEdgeSources_3;
                    				afterParserOrEnumRuleCall();
                    			

                    }


                    }
                    break;
                case 3 :
                    // InternalElkGraphJson.g:758:3: ( ruleKeyTargets otherlv_5= ':' this_ElkEdgeTargets_6= ruleElkEdgeTargets[$current] )
                    {
                    // InternalElkGraphJson.g:758:3: ( ruleKeyTargets otherlv_5= ':' this_ElkEdgeTargets_6= ruleElkEdgeTargets[$current] )
                    // InternalElkGraphJson.g:759:4: ruleKeyTargets otherlv_5= ':' this_ElkEdgeTargets_6= ruleElkEdgeTargets[$current]
                    {

                    				newCompositeNode(grammarAccess.getEdgeElementAccess().getKeyTargetsParserRuleCall_2_0());
                    			
                    pushFollow(FOLLOW_7);
                    ruleKeyTargets();

                    state._fsp--;


                    				afterParserOrEnumRuleCall();
                    			
                    otherlv_5=(Token)match(input,16,FOLLOW_8); 

                    				newLeafNode(otherlv_5, grammarAccess.getEdgeElementAccess().getColonKeyword_2_1());
                    			

                    				if (current==null) {
                    					current = createModelElement(grammarAccess.getEdgeElementRule());
                    				}
                    				newCompositeNode(grammarAccess.getEdgeElementAccess().getElkEdgeTargetsParserRuleCall_2_2());
                    			
                    pushFollow(FOLLOW_2);
                    this_ElkEdgeTargets_6=ruleElkEdgeTargets(current);

                    state._fsp--;


                    				current = this_ElkEdgeTargets_6;
                    				afterParserOrEnumRuleCall();
                    			

                    }


                    }
                    break;
                case 4 :
                    // InternalElkGraphJson.g:783:3: ( ruleKeyLabels otherlv_8= ':' this_ElkGraphElementLabels_9= ruleElkGraphElementLabels[$current] )
                    {
                    // InternalElkGraphJson.g:783:3: ( ruleKeyLabels otherlv_8= ':' this_ElkGraphElementLabels_9= ruleElkGraphElementLabels[$current] )
                    // InternalElkGraphJson.g:784:4: ruleKeyLabels otherlv_8= ':' this_ElkGraphElementLabels_9= ruleElkGraphElementLabels[$current]
                    {

                    				newCompositeNode(grammarAccess.getEdgeElementAccess().getKeyLabelsParserRuleCall_3_0());
                    			
                    pushFollow(FOLLOW_7);
                    ruleKeyLabels();

                    state._fsp--;


                    				afterParserOrEnumRuleCall();
                    			
                    otherlv_8=(Token)match(input,16,FOLLOW_8); 

                    				newLeafNode(otherlv_8, grammarAccess.getEdgeElementAccess().getColonKeyword_3_1());
                    			

                    				if (current==null) {
                    					current = createModelElement(grammarAccess.getEdgeElementRule());
                    				}
                    				newCompositeNode(grammarAccess.getEdgeElementAccess().getElkGraphElementLabelsParserRuleCall_3_2());
                    			
                    pushFollow(FOLLOW_2);
                    this_ElkGraphElementLabels_9=ruleElkGraphElementLabels(current);

                    state._fsp--;


                    				current = this_ElkGraphElementLabels_9;
                    				afterParserOrEnumRuleCall();
                    			

                    }


                    }
                    break;
                case 5 :
                    // InternalElkGraphJson.g:808:3: ( ruleKeyLayoutOptions otherlv_11= ':' this_ElkGraphElementProperties_12= ruleElkGraphElementProperties[$current] )
                    {
                    // InternalElkGraphJson.g:808:3: ( ruleKeyLayoutOptions otherlv_11= ':' this_ElkGraphElementProperties_12= ruleElkGraphElementProperties[$current] )
                    // InternalElkGraphJson.g:809:4: ruleKeyLayoutOptions otherlv_11= ':' this_ElkGraphElementProperties_12= ruleElkGraphElementProperties[$current]
                    {

                    				newCompositeNode(grammarAccess.getEdgeElementAccess().getKeyLayoutOptionsParserRuleCall_4_0());
                    			
                    pushFollow(FOLLOW_7);
                    ruleKeyLayoutOptions();

                    state._fsp--;


                    				afterParserOrEnumRuleCall();
                    			
                    otherlv_11=(Token)match(input,16,FOLLOW_9); 

                    				newLeafNode(otherlv_11, grammarAccess.getEdgeElementAccess().getColonKeyword_4_1());
                    			

                    				if (current==null) {
                    					current = createModelElement(grammarAccess.getEdgeElementRule());
                    				}
                    				newCompositeNode(grammarAccess.getEdgeElementAccess().getElkGraphElementPropertiesParserRuleCall_4_2());
                    			
                    pushFollow(FOLLOW_2);
                    this_ElkGraphElementProperties_12=ruleElkGraphElementProperties(current);

                    state._fsp--;


                    				current = this_ElkGraphElementProperties_12;
                    				afterParserOrEnumRuleCall();
                    			

                    }


                    }
                    break;
                case 6 :
                    // InternalElkGraphJson.g:833:3: ruleJsonMember
                    {

                    			newCompositeNode(grammarAccess.getEdgeElementAccess().getJsonMemberParserRuleCall_5());
                    		
                    pushFollow(FOLLOW_2);
                    ruleJsonMember();

                    state._fsp--;


                    			afterParserOrEnumRuleCall();
                    		

                    }
                    break;

            }


            }


            	leaveRule();

        }

            catch (RecognitionException re) {
                recover(input,re);
                appendSkippedTokens();
            }
        finally {
        }
        return current;
    }
    // $ANTLR end "ruleEdgeElement"


    // $ANTLR start "ruleElkEdgeSources"
    // InternalElkGraphJson.g:845:1: ruleElkEdgeSources[EObject in_current] returns [EObject current=in_current] : (otherlv_0= '[' ( ( (otherlv_1= RULE_STRING ) ) (otherlv_2= ',' ( (otherlv_3= RULE_STRING ) ) )* )? (otherlv_4= ',' )? otherlv_5= ']' ) ;
    public final EObject ruleElkEdgeSources(EObject in_current) throws RecognitionException {
        EObject current = in_current;

        Token otherlv_0=null;
        Token otherlv_1=null;
        Token otherlv_2=null;
        Token otherlv_3=null;
        Token otherlv_4=null;
        Token otherlv_5=null;


        	enterRule();

        try {
            // InternalElkGraphJson.g:851:2: ( (otherlv_0= '[' ( ( (otherlv_1= RULE_STRING ) ) (otherlv_2= ',' ( (otherlv_3= RULE_STRING ) ) )* )? (otherlv_4= ',' )? otherlv_5= ']' ) )
            // InternalElkGraphJson.g:852:2: (otherlv_0= '[' ( ( (otherlv_1= RULE_STRING ) ) (otherlv_2= ',' ( (otherlv_3= RULE_STRING ) ) )* )? (otherlv_4= ',' )? otherlv_5= ']' )
            {
            // InternalElkGraphJson.g:852:2: (otherlv_0= '[' ( ( (otherlv_1= RULE_STRING ) ) (otherlv_2= ',' ( (otherlv_3= RULE_STRING ) ) )* )? (otherlv_4= ',' )? otherlv_5= ']' )
            // InternalElkGraphJson.g:853:3: otherlv_0= '[' ( ( (otherlv_1= RULE_STRING ) ) (otherlv_2= ',' ( (otherlv_3= RULE_STRING ) ) )* )? (otherlv_4= ',' )? otherlv_5= ']'
            {
            otherlv_0=(Token)match(input,17,FOLLOW_13); 

            			newLeafNode(otherlv_0, grammarAccess.getElkEdgeSourcesAccess().getLeftSquareBracketKeyword_0());
            		
            // InternalElkGraphJson.g:857:3: ( ( (otherlv_1= RULE_STRING ) ) (otherlv_2= ',' ( (otherlv_3= RULE_STRING ) ) )* )?
            int alt15=2;
            int LA15_0 = input.LA(1);

            if ( (LA15_0==RULE_STRING) ) {
                alt15=1;
            }
            switch (alt15) {
                case 1 :
                    // InternalElkGraphJson.g:858:4: ( (otherlv_1= RULE_STRING ) ) (otherlv_2= ',' ( (otherlv_3= RULE_STRING ) ) )*
                    {
                    // InternalElkGraphJson.g:858:4: ( (otherlv_1= RULE_STRING ) )
                    // InternalElkGraphJson.g:859:5: (otherlv_1= RULE_STRING )
                    {
                    // InternalElkGraphJson.g:859:5: (otherlv_1= RULE_STRING )
                    // InternalElkGraphJson.g:860:6: otherlv_1= RULE_STRING
                    {

                    						if (current==null) {
                    							current = createModelElement(grammarAccess.getElkEdgeSourcesRule());
                    						}
                    					
                    otherlv_1=(Token)match(input,RULE_STRING,FOLLOW_14); 

                    						newLeafNode(otherlv_1, grammarAccess.getElkEdgeSourcesAccess().getSourcesElkConnectableShapeCrossReference_1_0_0());
                    					

                    }


                    }

                    // InternalElkGraphJson.g:871:4: (otherlv_2= ',' ( (otherlv_3= RULE_STRING ) ) )*
                    loop14:
                    do {
                        int alt14=2;
                        int LA14_0 = input.LA(1);

                        if ( (LA14_0==14) ) {
                            int LA14_1 = input.LA(2);

                            if ( (LA14_1==RULE_STRING) ) {
                                alt14=1;
                            }


                        }


                        switch (alt14) {
                    	case 1 :
                    	    // InternalElkGraphJson.g:872:5: otherlv_2= ',' ( (otherlv_3= RULE_STRING ) )
                    	    {
                    	    otherlv_2=(Token)match(input,14,FOLLOW_11); 

                    	    					newLeafNode(otherlv_2, grammarAccess.getElkEdgeSourcesAccess().getCommaKeyword_1_1_0());
                    	    				
                    	    // InternalElkGraphJson.g:876:5: ( (otherlv_3= RULE_STRING ) )
                    	    // InternalElkGraphJson.g:877:6: (otherlv_3= RULE_STRING )
                    	    {
                    	    // InternalElkGraphJson.g:877:6: (otherlv_3= RULE_STRING )
                    	    // InternalElkGraphJson.g:878:7: otherlv_3= RULE_STRING
                    	    {

                    	    							if (current==null) {
                    	    								current = createModelElement(grammarAccess.getElkEdgeSourcesRule());
                    	    							}
                    	    						
                    	    otherlv_3=(Token)match(input,RULE_STRING,FOLLOW_14); 

                    	    							newLeafNode(otherlv_3, grammarAccess.getElkEdgeSourcesAccess().getSourcesElkConnectableShapeCrossReference_1_1_1_0());
                    	    						

                    	    }


                    	    }


                    	    }
                    	    break;

                    	default :
                    	    break loop14;
                        }
                    } while (true);


                    }
                    break;

            }

            // InternalElkGraphJson.g:891:3: (otherlv_4= ',' )?
            int alt16=2;
            int LA16_0 = input.LA(1);

            if ( (LA16_0==14) ) {
                alt16=1;
            }
            switch (alt16) {
                case 1 :
                    // InternalElkGraphJson.g:892:4: otherlv_4= ','
                    {
                    otherlv_4=(Token)match(input,14,FOLLOW_15); 

                    				newLeafNode(otherlv_4, grammarAccess.getElkEdgeSourcesAccess().getCommaKeyword_2());
                    			

                    }
                    break;

            }

            otherlv_5=(Token)match(input,18,FOLLOW_2); 

            			newLeafNode(otherlv_5, grammarAccess.getElkEdgeSourcesAccess().getRightSquareBracketKeyword_3());
            		

            }


            }


            	leaveRule();

        }

            catch (RecognitionException re) {
                recover(input,re);
                appendSkippedTokens();
            }
        finally {
        }
        return current;
    }
    // $ANTLR end "ruleElkEdgeSources"


    // $ANTLR start "ruleElkEdgeTargets"
    // InternalElkGraphJson.g:906:1: ruleElkEdgeTargets[EObject in_current] returns [EObject current=in_current] : (otherlv_0= '[' ( ( (otherlv_1= RULE_STRING ) ) (otherlv_2= ',' ( (otherlv_3= RULE_STRING ) ) )* )? (otherlv_4= ',' )? otherlv_5= ']' ) ;
    public final EObject ruleElkEdgeTargets(EObject in_current) throws RecognitionException {
        EObject current = in_current;

        Token otherlv_0=null;
        Token otherlv_1=null;
        Token otherlv_2=null;
        Token otherlv_3=null;
        Token otherlv_4=null;
        Token otherlv_5=null;


        	enterRule();

        try {
            // InternalElkGraphJson.g:912:2: ( (otherlv_0= '[' ( ( (otherlv_1= RULE_STRING ) ) (otherlv_2= ',' ( (otherlv_3= RULE_STRING ) ) )* )? (otherlv_4= ',' )? otherlv_5= ']' ) )
            // InternalElkGraphJson.g:913:2: (otherlv_0= '[' ( ( (otherlv_1= RULE_STRING ) ) (otherlv_2= ',' ( (otherlv_3= RULE_STRING ) ) )* )? (otherlv_4= ',' )? otherlv_5= ']' )
            {
            // InternalElkGraphJson.g:913:2: (otherlv_0= '[' ( ( (otherlv_1= RULE_STRING ) ) (otherlv_2= ',' ( (otherlv_3= RULE_STRING ) ) )* )? (otherlv_4= ',' )? otherlv_5= ']' )
            // InternalElkGraphJson.g:914:3: otherlv_0= '[' ( ( (otherlv_1= RULE_STRING ) ) (otherlv_2= ',' ( (otherlv_3= RULE_STRING ) ) )* )? (otherlv_4= ',' )? otherlv_5= ']'
            {
            otherlv_0=(Token)match(input,17,FOLLOW_13); 

            			newLeafNode(otherlv_0, grammarAccess.getElkEdgeTargetsAccess().getLeftSquareBracketKeyword_0());
            		
            // InternalElkGraphJson.g:918:3: ( ( (otherlv_1= RULE_STRING ) ) (otherlv_2= ',' ( (otherlv_3= RULE_STRING ) ) )* )?
            int alt18=2;
            int LA18_0 = input.LA(1);

            if ( (LA18_0==RULE_STRING) ) {
                alt18=1;
            }
            switch (alt18) {
                case 1 :
                    // InternalElkGraphJson.g:919:4: ( (otherlv_1= RULE_STRING ) ) (otherlv_2= ',' ( (otherlv_3= RULE_STRING ) ) )*
                    {
                    // InternalElkGraphJson.g:919:4: ( (otherlv_1= RULE_STRING ) )
                    // InternalElkGraphJson.g:920:5: (otherlv_1= RULE_STRING )
                    {
                    // InternalElkGraphJson.g:920:5: (otherlv_1= RULE_STRING )
                    // InternalElkGraphJson.g:921:6: otherlv_1= RULE_STRING
                    {

                    						if (current==null) {
                    							current = createModelElement(grammarAccess.getElkEdgeTargetsRule());
                    						}
                    					
                    otherlv_1=(Token)match(input,RULE_STRING,FOLLOW_14); 

                    						newLeafNode(otherlv_1, grammarAccess.getElkEdgeTargetsAccess().getTargetsElkConnectableShapeCrossReference_1_0_0());
                    					

                    }


                    }

                    // InternalElkGraphJson.g:932:4: (otherlv_2= ',' ( (otherlv_3= RULE_STRING ) ) )*
                    loop17:
                    do {
                        int alt17=2;
                        int LA17_0 = input.LA(1);

                        if ( (LA17_0==14) ) {
                            int LA17_1 = input.LA(2);

                            if ( (LA17_1==RULE_STRING) ) {
                                alt17=1;
                            }


                        }


                        switch (alt17) {
                    	case 1 :
                    	    // InternalElkGraphJson.g:933:5: otherlv_2= ',' ( (otherlv_3= RULE_STRING ) )
                    	    {
                    	    otherlv_2=(Token)match(input,14,FOLLOW_11); 

                    	    					newLeafNode(otherlv_2, grammarAccess.getElkEdgeTargetsAccess().getCommaKeyword_1_1_0());
                    	    				
                    	    // InternalElkGraphJson.g:937:5: ( (otherlv_3= RULE_STRING ) )
                    	    // InternalElkGraphJson.g:938:6: (otherlv_3= RULE_STRING )
                    	    {
                    	    // InternalElkGraphJson.g:938:6: (otherlv_3= RULE_STRING )
                    	    // InternalElkGraphJson.g:939:7: otherlv_3= RULE_STRING
                    	    {

                    	    							if (current==null) {
                    	    								current = createModelElement(grammarAccess.getElkEdgeTargetsRule());
                    	    							}
                    	    						
                    	    otherlv_3=(Token)match(input,RULE_STRING,FOLLOW_14); 

                    	    							newLeafNode(otherlv_3, grammarAccess.getElkEdgeTargetsAccess().getTargetsElkConnectableShapeCrossReference_1_1_1_0());
                    	    						

                    	    }


                    	    }


                    	    }
                    	    break;

                    	default :
                    	    break loop17;
                        }
                    } while (true);


                    }
                    break;

            }

            // InternalElkGraphJson.g:952:3: (otherlv_4= ',' )?
            int alt19=2;
            int LA19_0 = input.LA(1);

            if ( (LA19_0==14) ) {
                alt19=1;
            }
            switch (alt19) {
                case 1 :
                    // InternalElkGraphJson.g:953:4: otherlv_4= ','
                    {
                    otherlv_4=(Token)match(input,14,FOLLOW_15); 

                    				newLeafNode(otherlv_4, grammarAccess.getElkEdgeTargetsAccess().getCommaKeyword_2());
                    			

                    }
                    break;

            }

            otherlv_5=(Token)match(input,18,FOLLOW_2); 

            			newLeafNode(otherlv_5, grammarAccess.getElkEdgeTargetsAccess().getRightSquareBracketKeyword_3());
            		

            }


            }


            	leaveRule();

        }

            catch (RecognitionException re) {
                recover(input,re);
                appendSkippedTokens();
            }
        finally {
        }
        return current;
    }
    // $ANTLR end "ruleElkEdgeTargets"


    // $ANTLR start "ruleElkId"
    // InternalElkGraphJson.g:967:1: ruleElkId[EObject in_current] returns [EObject current=in_current] : ( ruleKeyId otherlv_1= ':' ( (lv_identifier_2_0= RULE_STRING ) ) ) ;
    public final EObject ruleElkId(EObject in_current) throws RecognitionException {
        EObject current = in_current;

        Token otherlv_1=null;
        Token lv_identifier_2_0=null;


        	enterRule();

        try {
            // InternalElkGraphJson.g:973:2: ( ( ruleKeyId otherlv_1= ':' ( (lv_identifier_2_0= RULE_STRING ) ) ) )
            // InternalElkGraphJson.g:974:2: ( ruleKeyId otherlv_1= ':' ( (lv_identifier_2_0= RULE_STRING ) ) )
            {
            // InternalElkGraphJson.g:974:2: ( ruleKeyId otherlv_1= ':' ( (lv_identifier_2_0= RULE_STRING ) ) )
            // InternalElkGraphJson.g:975:3: ruleKeyId otherlv_1= ':' ( (lv_identifier_2_0= RULE_STRING ) )
            {

            			newCompositeNode(grammarAccess.getElkIdAccess().getKeyIdParserRuleCall_0());
            		
            pushFollow(FOLLOW_7);
            ruleKeyId();

            state._fsp--;


            			afterParserOrEnumRuleCall();
            		
            otherlv_1=(Token)match(input,16,FOLLOW_11); 

            			newLeafNode(otherlv_1, grammarAccess.getElkIdAccess().getColonKeyword_1());
            		
            // InternalElkGraphJson.g:986:3: ( (lv_identifier_2_0= RULE_STRING ) )
            // InternalElkGraphJson.g:987:4: (lv_identifier_2_0= RULE_STRING )
            {
            // InternalElkGraphJson.g:987:4: (lv_identifier_2_0= RULE_STRING )
            // InternalElkGraphJson.g:988:5: lv_identifier_2_0= RULE_STRING
            {
            lv_identifier_2_0=(Token)match(input,RULE_STRING,FOLLOW_2); 

            					newLeafNode(lv_identifier_2_0, grammarAccess.getElkIdAccess().getIdentifierSTRINGTerminalRuleCall_2_0());
            				

            					if (current==null) {
            						current = createModelElement(grammarAccess.getElkIdRule());
            					}
            					setWithLastConsumed(
            						current,
            						"identifier",
            						lv_identifier_2_0,
            						"org.eclipse.xtext.common.Terminals.STRING");
            				

            }


            }


            }


            }


            	leaveRule();

        }

            catch (RecognitionException re) {
                recover(input,re);
                appendSkippedTokens();
            }
        finally {
        }
        return current;
    }
    // $ANTLR end "ruleElkId"


    // $ANTLR start "ruleElkNodeChildren"
    // InternalElkGraphJson.g:1009:1: ruleElkNodeChildren[EObject in_current] returns [EObject current=in_current] : (otherlv_0= '[' ( ( (lv_children_1_0= ruleElkNode ) ) (otherlv_2= ',' ( (lv_children_3_0= ruleElkNode ) ) )* )? (otherlv_4= ',' )? otherlv_5= ']' ) ;
    public final EObject ruleElkNodeChildren(EObject in_current) throws RecognitionException {
        EObject current = in_current;

        Token otherlv_0=null;
        Token otherlv_2=null;
        Token otherlv_4=null;
        Token otherlv_5=null;
        EObject lv_children_1_0 = null;

        EObject lv_children_3_0 = null;



        	enterRule();

        try {
            // InternalElkGraphJson.g:1015:2: ( (otherlv_0= '[' ( ( (lv_children_1_0= ruleElkNode ) ) (otherlv_2= ',' ( (lv_children_3_0= ruleElkNode ) ) )* )? (otherlv_4= ',' )? otherlv_5= ']' ) )
            // InternalElkGraphJson.g:1016:2: (otherlv_0= '[' ( ( (lv_children_1_0= ruleElkNode ) ) (otherlv_2= ',' ( (lv_children_3_0= ruleElkNode ) ) )* )? (otherlv_4= ',' )? otherlv_5= ']' )
            {
            // InternalElkGraphJson.g:1016:2: (otherlv_0= '[' ( ( (lv_children_1_0= ruleElkNode ) ) (otherlv_2= ',' ( (lv_children_3_0= ruleElkNode ) ) )* )? (otherlv_4= ',' )? otherlv_5= ']' )
            // InternalElkGraphJson.g:1017:3: otherlv_0= '[' ( ( (lv_children_1_0= ruleElkNode ) ) (otherlv_2= ',' ( (lv_children_3_0= ruleElkNode ) ) )* )? (otherlv_4= ',' )? otherlv_5= ']'
            {
            otherlv_0=(Token)match(input,17,FOLLOW_16); 

            			newLeafNode(otherlv_0, grammarAccess.getElkNodeChildrenAccess().getLeftSquareBracketKeyword_0());
            		
            // InternalElkGraphJson.g:1021:3: ( ( (lv_children_1_0= ruleElkNode ) ) (otherlv_2= ',' ( (lv_children_3_0= ruleElkNode ) ) )* )?
            int alt21=2;
            int LA21_0 = input.LA(1);

            if ( (LA21_0==13) ) {
                alt21=1;
            }
            switch (alt21) {
                case 1 :
                    // InternalElkGraphJson.g:1022:4: ( (lv_children_1_0= ruleElkNode ) ) (otherlv_2= ',' ( (lv_children_3_0= ruleElkNode ) ) )*
                    {
                    // InternalElkGraphJson.g:1022:4: ( (lv_children_1_0= ruleElkNode ) )
                    // InternalElkGraphJson.g:1023:5: (lv_children_1_0= ruleElkNode )
                    {
                    // InternalElkGraphJson.g:1023:5: (lv_children_1_0= ruleElkNode )
                    // InternalElkGraphJson.g:1024:6: lv_children_1_0= ruleElkNode
                    {

                    						newCompositeNode(grammarAccess.getElkNodeChildrenAccess().getChildrenElkNodeParserRuleCall_1_0_0());
                    					
                    pushFollow(FOLLOW_14);
                    lv_children_1_0=ruleElkNode();

                    state._fsp--;


                    						if (current==null) {
                    							current = createModelElementForParent(grammarAccess.getElkNodeChildrenRule());
                    						}
                    						add(
                    							current,
                    							"children",
                    							lv_children_1_0,
                    							"org.eclipse.elk.graph.json.text.ElkGraphJson.ElkNode");
                    						afterParserOrEnumRuleCall();
                    					

                    }


                    }

                    // InternalElkGraphJson.g:1041:4: (otherlv_2= ',' ( (lv_children_3_0= ruleElkNode ) ) )*
                    loop20:
                    do {
                        int alt20=2;
                        int LA20_0 = input.LA(1);

                        if ( (LA20_0==14) ) {
                            int LA20_1 = input.LA(2);

                            if ( (LA20_1==13) ) {
                                alt20=1;
                            }


                        }


                        switch (alt20) {
                    	case 1 :
                    	    // InternalElkGraphJson.g:1042:5: otherlv_2= ',' ( (lv_children_3_0= ruleElkNode ) )
                    	    {
                    	    otherlv_2=(Token)match(input,14,FOLLOW_9); 

                    	    					newLeafNode(otherlv_2, grammarAccess.getElkNodeChildrenAccess().getCommaKeyword_1_1_0());
                    	    				
                    	    // InternalElkGraphJson.g:1046:5: ( (lv_children_3_0= ruleElkNode ) )
                    	    // InternalElkGraphJson.g:1047:6: (lv_children_3_0= ruleElkNode )
                    	    {
                    	    // InternalElkGraphJson.g:1047:6: (lv_children_3_0= ruleElkNode )
                    	    // InternalElkGraphJson.g:1048:7: lv_children_3_0= ruleElkNode
                    	    {

                    	    							newCompositeNode(grammarAccess.getElkNodeChildrenAccess().getChildrenElkNodeParserRuleCall_1_1_1_0());
                    	    						
                    	    pushFollow(FOLLOW_14);
                    	    lv_children_3_0=ruleElkNode();

                    	    state._fsp--;


                    	    							if (current==null) {
                    	    								current = createModelElementForParent(grammarAccess.getElkNodeChildrenRule());
                    	    							}
                    	    							add(
                    	    								current,
                    	    								"children",
                    	    								lv_children_3_0,
                    	    								"org.eclipse.elk.graph.json.text.ElkGraphJson.ElkNode");
                    	    							afterParserOrEnumRuleCall();
                    	    						

                    	    }


                    	    }


                    	    }
                    	    break;

                    	default :
                    	    break loop20;
                        }
                    } while (true);


                    }
                    break;

            }

            // InternalElkGraphJson.g:1067:3: (otherlv_4= ',' )?
            int alt22=2;
            int LA22_0 = input.LA(1);

            if ( (LA22_0==14) ) {
                alt22=1;
            }
            switch (alt22) {
                case 1 :
                    // InternalElkGraphJson.g:1068:4: otherlv_4= ','
                    {
                    otherlv_4=(Token)match(input,14,FOLLOW_15); 

                    				newLeafNode(otherlv_4, grammarAccess.getElkNodeChildrenAccess().getCommaKeyword_2());
                    			

                    }
                    break;

            }

            otherlv_5=(Token)match(input,18,FOLLOW_2); 

            			newLeafNode(otherlv_5, grammarAccess.getElkNodeChildrenAccess().getRightSquareBracketKeyword_3());
            		

            }


            }


            	leaveRule();

        }

            catch (RecognitionException re) {
                recover(input,re);
                appendSkippedTokens();
            }
        finally {
        }
        return current;
    }
    // $ANTLR end "ruleElkNodeChildren"


    // $ANTLR start "ruleElkNodePorts"
    // InternalElkGraphJson.g:1082:1: ruleElkNodePorts[EObject in_current] returns [EObject current=in_current] : (otherlv_0= '[' ( ( (lv_ports_1_0= ruleElkPort ) ) (otherlv_2= ',' ( (lv_ports_3_0= ruleElkPort ) ) )* )? (otherlv_4= ',' )? otherlv_5= ']' ) ;
    public final EObject ruleElkNodePorts(EObject in_current) throws RecognitionException {
        EObject current = in_current;

        Token otherlv_0=null;
        Token otherlv_2=null;
        Token otherlv_4=null;
        Token otherlv_5=null;
        EObject lv_ports_1_0 = null;

        EObject lv_ports_3_0 = null;



        	enterRule();

        try {
            // InternalElkGraphJson.g:1088:2: ( (otherlv_0= '[' ( ( (lv_ports_1_0= ruleElkPort ) ) (otherlv_2= ',' ( (lv_ports_3_0= ruleElkPort ) ) )* )? (otherlv_4= ',' )? otherlv_5= ']' ) )
            // InternalElkGraphJson.g:1089:2: (otherlv_0= '[' ( ( (lv_ports_1_0= ruleElkPort ) ) (otherlv_2= ',' ( (lv_ports_3_0= ruleElkPort ) ) )* )? (otherlv_4= ',' )? otherlv_5= ']' )
            {
            // InternalElkGraphJson.g:1089:2: (otherlv_0= '[' ( ( (lv_ports_1_0= ruleElkPort ) ) (otherlv_2= ',' ( (lv_ports_3_0= ruleElkPort ) ) )* )? (otherlv_4= ',' )? otherlv_5= ']' )
            // InternalElkGraphJson.g:1090:3: otherlv_0= '[' ( ( (lv_ports_1_0= ruleElkPort ) ) (otherlv_2= ',' ( (lv_ports_3_0= ruleElkPort ) ) )* )? (otherlv_4= ',' )? otherlv_5= ']'
            {
            otherlv_0=(Token)match(input,17,FOLLOW_16); 

            			newLeafNode(otherlv_0, grammarAccess.getElkNodePortsAccess().getLeftSquareBracketKeyword_0());
            		
            // InternalElkGraphJson.g:1094:3: ( ( (lv_ports_1_0= ruleElkPort ) ) (otherlv_2= ',' ( (lv_ports_3_0= ruleElkPort ) ) )* )?
            int alt24=2;
            int LA24_0 = input.LA(1);

            if ( (LA24_0==13) ) {
                alt24=1;
            }
            switch (alt24) {
                case 1 :
                    // InternalElkGraphJson.g:1095:4: ( (lv_ports_1_0= ruleElkPort ) ) (otherlv_2= ',' ( (lv_ports_3_0= ruleElkPort ) ) )*
                    {
                    // InternalElkGraphJson.g:1095:4: ( (lv_ports_1_0= ruleElkPort ) )
                    // InternalElkGraphJson.g:1096:5: (lv_ports_1_0= ruleElkPort )
                    {
                    // InternalElkGraphJson.g:1096:5: (lv_ports_1_0= ruleElkPort )
                    // InternalElkGraphJson.g:1097:6: lv_ports_1_0= ruleElkPort
                    {

                    						newCompositeNode(grammarAccess.getElkNodePortsAccess().getPortsElkPortParserRuleCall_1_0_0());
                    					
                    pushFollow(FOLLOW_14);
                    lv_ports_1_0=ruleElkPort();

                    state._fsp--;


                    						if (current==null) {
                    							current = createModelElementForParent(grammarAccess.getElkNodePortsRule());
                    						}
                    						add(
                    							current,
                    							"ports",
                    							lv_ports_1_0,
                    							"org.eclipse.elk.graph.json.text.ElkGraphJson.ElkPort");
                    						afterParserOrEnumRuleCall();
                    					

                    }


                    }

                    // InternalElkGraphJson.g:1114:4: (otherlv_2= ',' ( (lv_ports_3_0= ruleElkPort ) ) )*
                    loop23:
                    do {
                        int alt23=2;
                        int LA23_0 = input.LA(1);

                        if ( (LA23_0==14) ) {
                            int LA23_1 = input.LA(2);

                            if ( (LA23_1==13) ) {
                                alt23=1;
                            }


                        }


                        switch (alt23) {
                    	case 1 :
                    	    // InternalElkGraphJson.g:1115:5: otherlv_2= ',' ( (lv_ports_3_0= ruleElkPort ) )
                    	    {
                    	    otherlv_2=(Token)match(input,14,FOLLOW_9); 

                    	    					newLeafNode(otherlv_2, grammarAccess.getElkNodePortsAccess().getCommaKeyword_1_1_0());
                    	    				
                    	    // InternalElkGraphJson.g:1119:5: ( (lv_ports_3_0= ruleElkPort ) )
                    	    // InternalElkGraphJson.g:1120:6: (lv_ports_3_0= ruleElkPort )
                    	    {
                    	    // InternalElkGraphJson.g:1120:6: (lv_ports_3_0= ruleElkPort )
                    	    // InternalElkGraphJson.g:1121:7: lv_ports_3_0= ruleElkPort
                    	    {

                    	    							newCompositeNode(grammarAccess.getElkNodePortsAccess().getPortsElkPortParserRuleCall_1_1_1_0());
                    	    						
                    	    pushFollow(FOLLOW_14);
                    	    lv_ports_3_0=ruleElkPort();

                    	    state._fsp--;


                    	    							if (current==null) {
                    	    								current = createModelElementForParent(grammarAccess.getElkNodePortsRule());
                    	    							}
                    	    							add(
                    	    								current,
                    	    								"ports",
                    	    								lv_ports_3_0,
                    	    								"org.eclipse.elk.graph.json.text.ElkGraphJson.ElkPort");
                    	    							afterParserOrEnumRuleCall();
                    	    						

                    	    }


                    	    }


                    	    }
                    	    break;

                    	default :
                    	    break loop23;
                        }
                    } while (true);


                    }
                    break;

            }

            // InternalElkGraphJson.g:1140:3: (otherlv_4= ',' )?
            int alt25=2;
            int LA25_0 = input.LA(1);

            if ( (LA25_0==14) ) {
                alt25=1;
            }
            switch (alt25) {
                case 1 :
                    // InternalElkGraphJson.g:1141:4: otherlv_4= ','
                    {
                    otherlv_4=(Token)match(input,14,FOLLOW_15); 

                    				newLeafNode(otherlv_4, grammarAccess.getElkNodePortsAccess().getCommaKeyword_2());
                    			

                    }
                    break;

            }

            otherlv_5=(Token)match(input,18,FOLLOW_2); 

            			newLeafNode(otherlv_5, grammarAccess.getElkNodePortsAccess().getRightSquareBracketKeyword_3());
            		

            }


            }


            	leaveRule();

        }

            catch (RecognitionException re) {
                recover(input,re);
                appendSkippedTokens();
            }
        finally {
        }
        return current;
    }
    // $ANTLR end "ruleElkNodePorts"


    // $ANTLR start "ruleElkNodeEdges"
    // InternalElkGraphJson.g:1155:1: ruleElkNodeEdges[EObject in_current] returns [EObject current=in_current] : (otherlv_0= '[' ( ( (lv_containedEdges_1_0= ruleElkEdge ) ) (otherlv_2= ',' ( (lv_containedEdges_3_0= ruleElkEdge ) ) )* )? (otherlv_4= ',' )? otherlv_5= ']' ) ;
    public final EObject ruleElkNodeEdges(EObject in_current) throws RecognitionException {
        EObject current = in_current;

        Token otherlv_0=null;
        Token otherlv_2=null;
        Token otherlv_4=null;
        Token otherlv_5=null;
        EObject lv_containedEdges_1_0 = null;

        EObject lv_containedEdges_3_0 = null;



        	enterRule();

        try {
            // InternalElkGraphJson.g:1161:2: ( (otherlv_0= '[' ( ( (lv_containedEdges_1_0= ruleElkEdge ) ) (otherlv_2= ',' ( (lv_containedEdges_3_0= ruleElkEdge ) ) )* )? (otherlv_4= ',' )? otherlv_5= ']' ) )
            // InternalElkGraphJson.g:1162:2: (otherlv_0= '[' ( ( (lv_containedEdges_1_0= ruleElkEdge ) ) (otherlv_2= ',' ( (lv_containedEdges_3_0= ruleElkEdge ) ) )* )? (otherlv_4= ',' )? otherlv_5= ']' )
            {
            // InternalElkGraphJson.g:1162:2: (otherlv_0= '[' ( ( (lv_containedEdges_1_0= ruleElkEdge ) ) (otherlv_2= ',' ( (lv_containedEdges_3_0= ruleElkEdge ) ) )* )? (otherlv_4= ',' )? otherlv_5= ']' )
            // InternalElkGraphJson.g:1163:3: otherlv_0= '[' ( ( (lv_containedEdges_1_0= ruleElkEdge ) ) (otherlv_2= ',' ( (lv_containedEdges_3_0= ruleElkEdge ) ) )* )? (otherlv_4= ',' )? otherlv_5= ']'
            {
            otherlv_0=(Token)match(input,17,FOLLOW_16); 

            			newLeafNode(otherlv_0, grammarAccess.getElkNodeEdgesAccess().getLeftSquareBracketKeyword_0());
            		
            // InternalElkGraphJson.g:1167:3: ( ( (lv_containedEdges_1_0= ruleElkEdge ) ) (otherlv_2= ',' ( (lv_containedEdges_3_0= ruleElkEdge ) ) )* )?
            int alt27=2;
            int LA27_0 = input.LA(1);

            if ( (LA27_0==13) ) {
                alt27=1;
            }
            switch (alt27) {
                case 1 :
                    // InternalElkGraphJson.g:1168:4: ( (lv_containedEdges_1_0= ruleElkEdge ) ) (otherlv_2= ',' ( (lv_containedEdges_3_0= ruleElkEdge ) ) )*
                    {
                    // InternalElkGraphJson.g:1168:4: ( (lv_containedEdges_1_0= ruleElkEdge ) )
                    // InternalElkGraphJson.g:1169:5: (lv_containedEdges_1_0= ruleElkEdge )
                    {
                    // InternalElkGraphJson.g:1169:5: (lv_containedEdges_1_0= ruleElkEdge )
                    // InternalElkGraphJson.g:1170:6: lv_containedEdges_1_0= ruleElkEdge
                    {

                    						newCompositeNode(grammarAccess.getElkNodeEdgesAccess().getContainedEdgesElkEdgeParserRuleCall_1_0_0());
                    					
                    pushFollow(FOLLOW_14);
                    lv_containedEdges_1_0=ruleElkEdge();

                    state._fsp--;


                    						if (current==null) {
                    							current = createModelElementForParent(grammarAccess.getElkNodeEdgesRule());
                    						}
                    						add(
                    							current,
                    							"containedEdges",
                    							lv_containedEdges_1_0,
                    							"org.eclipse.elk.graph.json.text.ElkGraphJson.ElkEdge");
                    						afterParserOrEnumRuleCall();
                    					

                    }


                    }

                    // InternalElkGraphJson.g:1187:4: (otherlv_2= ',' ( (lv_containedEdges_3_0= ruleElkEdge ) ) )*
                    loop26:
                    do {
                        int alt26=2;
                        int LA26_0 = input.LA(1);

                        if ( (LA26_0==14) ) {
                            int LA26_1 = input.LA(2);

                            if ( (LA26_1==13) ) {
                                alt26=1;
                            }


                        }


                        switch (alt26) {
                    	case 1 :
                    	    // InternalElkGraphJson.g:1188:5: otherlv_2= ',' ( (lv_containedEdges_3_0= ruleElkEdge ) )
                    	    {
                    	    otherlv_2=(Token)match(input,14,FOLLOW_9); 

                    	    					newLeafNode(otherlv_2, grammarAccess.getElkNodeEdgesAccess().getCommaKeyword_1_1_0());
                    	    				
                    	    // InternalElkGraphJson.g:1192:5: ( (lv_containedEdges_3_0= ruleElkEdge ) )
                    	    // InternalElkGraphJson.g:1193:6: (lv_containedEdges_3_0= ruleElkEdge )
                    	    {
                    	    // InternalElkGraphJson.g:1193:6: (lv_containedEdges_3_0= ruleElkEdge )
                    	    // InternalElkGraphJson.g:1194:7: lv_containedEdges_3_0= ruleElkEdge
                    	    {

                    	    							newCompositeNode(grammarAccess.getElkNodeEdgesAccess().getContainedEdgesElkEdgeParserRuleCall_1_1_1_0());
                    	    						
                    	    pushFollow(FOLLOW_14);
                    	    lv_containedEdges_3_0=ruleElkEdge();

                    	    state._fsp--;


                    	    							if (current==null) {
                    	    								current = createModelElementForParent(grammarAccess.getElkNodeEdgesRule());
                    	    							}
                    	    							add(
                    	    								current,
                    	    								"containedEdges",
                    	    								lv_containedEdges_3_0,
                    	    								"org.eclipse.elk.graph.json.text.ElkGraphJson.ElkEdge");
                    	    							afterParserOrEnumRuleCall();
                    	    						

                    	    }


                    	    }


                    	    }
                    	    break;

                    	default :
                    	    break loop26;
                        }
                    } while (true);


                    }
                    break;

            }

            // InternalElkGraphJson.g:1213:3: (otherlv_4= ',' )?
            int alt28=2;
            int LA28_0 = input.LA(1);

            if ( (LA28_0==14) ) {
                alt28=1;
            }
            switch (alt28) {
                case 1 :
                    // InternalElkGraphJson.g:1214:4: otherlv_4= ','
                    {
                    otherlv_4=(Token)match(input,14,FOLLOW_15); 

                    				newLeafNode(otherlv_4, grammarAccess.getElkNodeEdgesAccess().getCommaKeyword_2());
                    			

                    }
                    break;

            }

            otherlv_5=(Token)match(input,18,FOLLOW_2); 

            			newLeafNode(otherlv_5, grammarAccess.getElkNodeEdgesAccess().getRightSquareBracketKeyword_3());
            		

            }


            }


            	leaveRule();

        }

            catch (RecognitionException re) {
                recover(input,re);
                appendSkippedTokens();
            }
        finally {
        }
        return current;
    }
    // $ANTLR end "ruleElkNodeEdges"


    // $ANTLR start "ruleElkGraphElementLabels"
    // InternalElkGraphJson.g:1228:1: ruleElkGraphElementLabels[EObject in_current] returns [EObject current=in_current] : (otherlv_0= '[' ( ( (lv_labels_1_0= ruleElkLabel ) ) (otherlv_2= ',' ( (lv_labels_3_0= ruleElkLabel ) ) )* )? (otherlv_4= ',' )? otherlv_5= ']' ) ;
    public final EObject ruleElkGraphElementLabels(EObject in_current) throws RecognitionException {
        EObject current = in_current;

        Token otherlv_0=null;
        Token otherlv_2=null;
        Token otherlv_4=null;
        Token otherlv_5=null;
        EObject lv_labels_1_0 = null;

        EObject lv_labels_3_0 = null;



        	enterRule();

        try {
            // InternalElkGraphJson.g:1234:2: ( (otherlv_0= '[' ( ( (lv_labels_1_0= ruleElkLabel ) ) (otherlv_2= ',' ( (lv_labels_3_0= ruleElkLabel ) ) )* )? (otherlv_4= ',' )? otherlv_5= ']' ) )
            // InternalElkGraphJson.g:1235:2: (otherlv_0= '[' ( ( (lv_labels_1_0= ruleElkLabel ) ) (otherlv_2= ',' ( (lv_labels_3_0= ruleElkLabel ) ) )* )? (otherlv_4= ',' )? otherlv_5= ']' )
            {
            // InternalElkGraphJson.g:1235:2: (otherlv_0= '[' ( ( (lv_labels_1_0= ruleElkLabel ) ) (otherlv_2= ',' ( (lv_labels_3_0= ruleElkLabel ) ) )* )? (otherlv_4= ',' )? otherlv_5= ']' )
            // InternalElkGraphJson.g:1236:3: otherlv_0= '[' ( ( (lv_labels_1_0= ruleElkLabel ) ) (otherlv_2= ',' ( (lv_labels_3_0= ruleElkLabel ) ) )* )? (otherlv_4= ',' )? otherlv_5= ']'
            {
            otherlv_0=(Token)match(input,17,FOLLOW_16); 

            			newLeafNode(otherlv_0, grammarAccess.getElkGraphElementLabelsAccess().getLeftSquareBracketKeyword_0());
            		
            // InternalElkGraphJson.g:1240:3: ( ( (lv_labels_1_0= ruleElkLabel ) ) (otherlv_2= ',' ( (lv_labels_3_0= ruleElkLabel ) ) )* )?
            int alt30=2;
            int LA30_0 = input.LA(1);

            if ( (LA30_0==13) ) {
                alt30=1;
            }
            switch (alt30) {
                case 1 :
                    // InternalElkGraphJson.g:1241:4: ( (lv_labels_1_0= ruleElkLabel ) ) (otherlv_2= ',' ( (lv_labels_3_0= ruleElkLabel ) ) )*
                    {
                    // InternalElkGraphJson.g:1241:4: ( (lv_labels_1_0= ruleElkLabel ) )
                    // InternalElkGraphJson.g:1242:5: (lv_labels_1_0= ruleElkLabel )
                    {
                    // InternalElkGraphJson.g:1242:5: (lv_labels_1_0= ruleElkLabel )
                    // InternalElkGraphJson.g:1243:6: lv_labels_1_0= ruleElkLabel
                    {

                    						newCompositeNode(grammarAccess.getElkGraphElementLabelsAccess().getLabelsElkLabelParserRuleCall_1_0_0());
                    					
                    pushFollow(FOLLOW_14);
                    lv_labels_1_0=ruleElkLabel();

                    state._fsp--;


                    						if (current==null) {
                    							current = createModelElementForParent(grammarAccess.getElkGraphElementLabelsRule());
                    						}
                    						add(
                    							current,
                    							"labels",
                    							lv_labels_1_0,
                    							"org.eclipse.elk.graph.json.text.ElkGraphJson.ElkLabel");
                    						afterParserOrEnumRuleCall();
                    					

                    }


                    }

                    // InternalElkGraphJson.g:1260:4: (otherlv_2= ',' ( (lv_labels_3_0= ruleElkLabel ) ) )*
                    loop29:
                    do {
                        int alt29=2;
                        int LA29_0 = input.LA(1);

                        if ( (LA29_0==14) ) {
                            int LA29_1 = input.LA(2);

                            if ( (LA29_1==13) ) {
                                alt29=1;
                            }


                        }


                        switch (alt29) {
                    	case 1 :
                    	    // InternalElkGraphJson.g:1261:5: otherlv_2= ',' ( (lv_labels_3_0= ruleElkLabel ) )
                    	    {
                    	    otherlv_2=(Token)match(input,14,FOLLOW_9); 

                    	    					newLeafNode(otherlv_2, grammarAccess.getElkGraphElementLabelsAccess().getCommaKeyword_1_1_0());
                    	    				
                    	    // InternalElkGraphJson.g:1265:5: ( (lv_labels_3_0= ruleElkLabel ) )
                    	    // InternalElkGraphJson.g:1266:6: (lv_labels_3_0= ruleElkLabel )
                    	    {
                    	    // InternalElkGraphJson.g:1266:6: (lv_labels_3_0= ruleElkLabel )
                    	    // InternalElkGraphJson.g:1267:7: lv_labels_3_0= ruleElkLabel
                    	    {

                    	    							newCompositeNode(grammarAccess.getElkGraphElementLabelsAccess().getLabelsElkLabelParserRuleCall_1_1_1_0());
                    	    						
                    	    pushFollow(FOLLOW_14);
                    	    lv_labels_3_0=ruleElkLabel();

                    	    state._fsp--;


                    	    							if (current==null) {
                    	    								current = createModelElementForParent(grammarAccess.getElkGraphElementLabelsRule());
                    	    							}
                    	    							add(
                    	    								current,
                    	    								"labels",
                    	    								lv_labels_3_0,
                    	    								"org.eclipse.elk.graph.json.text.ElkGraphJson.ElkLabel");
                    	    							afterParserOrEnumRuleCall();
                    	    						

                    	    }


                    	    }


                    	    }
                    	    break;

                    	default :
                    	    break loop29;
                        }
                    } while (true);


                    }
                    break;

            }

            // InternalElkGraphJson.g:1286:3: (otherlv_4= ',' )?
            int alt31=2;
            int LA31_0 = input.LA(1);

            if ( (LA31_0==14) ) {
                alt31=1;
            }
            switch (alt31) {
                case 1 :
                    // InternalElkGraphJson.g:1287:4: otherlv_4= ','
                    {
                    otherlv_4=(Token)match(input,14,FOLLOW_15); 

                    				newLeafNode(otherlv_4, grammarAccess.getElkGraphElementLabelsAccess().getCommaKeyword_2());
                    			

                    }
                    break;

            }

            otherlv_5=(Token)match(input,18,FOLLOW_2); 

            			newLeafNode(otherlv_5, grammarAccess.getElkGraphElementLabelsAccess().getRightSquareBracketKeyword_3());
            		

            }


            }


            	leaveRule();

        }

            catch (RecognitionException re) {
                recover(input,re);
                appendSkippedTokens();
            }
        finally {
        }
        return current;
    }
    // $ANTLR end "ruleElkGraphElementLabels"


    // $ANTLR start "ruleElkGraphElementProperties"
    // InternalElkGraphJson.g:1301:1: ruleElkGraphElementProperties[EObject in_current] returns [EObject current=in_current] : (otherlv_0= '{' ( ( (lv_properties_1_0= ruleProperty ) ) (otherlv_2= ',' ( (lv_properties_3_0= ruleProperty ) ) )* )? (otherlv_4= ',' )? otherlv_5= '}' ) ;
    public final EObject ruleElkGraphElementProperties(EObject in_current) throws RecognitionException {
        EObject current = in_current;

        Token otherlv_0=null;
        Token otherlv_2=null;
        Token otherlv_4=null;
        Token otherlv_5=null;
        EObject lv_properties_1_0 = null;

        EObject lv_properties_3_0 = null;



        	enterRule();

        try {
            // InternalElkGraphJson.g:1307:2: ( (otherlv_0= '{' ( ( (lv_properties_1_0= ruleProperty ) ) (otherlv_2= ',' ( (lv_properties_3_0= ruleProperty ) ) )* )? (otherlv_4= ',' )? otherlv_5= '}' ) )
            // InternalElkGraphJson.g:1308:2: (otherlv_0= '{' ( ( (lv_properties_1_0= ruleProperty ) ) (otherlv_2= ',' ( (lv_properties_3_0= ruleProperty ) ) )* )? (otherlv_4= ',' )? otherlv_5= '}' )
            {
            // InternalElkGraphJson.g:1308:2: (otherlv_0= '{' ( ( (lv_properties_1_0= ruleProperty ) ) (otherlv_2= ',' ( (lv_properties_3_0= ruleProperty ) ) )* )? (otherlv_4= ',' )? otherlv_5= '}' )
            // InternalElkGraphJson.g:1309:3: otherlv_0= '{' ( ( (lv_properties_1_0= ruleProperty ) ) (otherlv_2= ',' ( (lv_properties_3_0= ruleProperty ) ) )* )? (otherlv_4= ',' )? otherlv_5= '}'
            {
            otherlv_0=(Token)match(input,13,FOLLOW_17); 

            			newLeafNode(otherlv_0, grammarAccess.getElkGraphElementPropertiesAccess().getLeftCurlyBracketKeyword_0());
            		
            // InternalElkGraphJson.g:1313:3: ( ( (lv_properties_1_0= ruleProperty ) ) (otherlv_2= ',' ( (lv_properties_3_0= ruleProperty ) ) )* )?
            int alt33=2;
            int LA33_0 = input.LA(1);

            if ( ((LA33_0>=RULE_STRING && LA33_0<=RULE_ID)) ) {
                alt33=1;
            }
            switch (alt33) {
                case 1 :
                    // InternalElkGraphJson.g:1314:4: ( (lv_properties_1_0= ruleProperty ) ) (otherlv_2= ',' ( (lv_properties_3_0= ruleProperty ) ) )*
                    {
                    // InternalElkGraphJson.g:1314:4: ( (lv_properties_1_0= ruleProperty ) )
                    // InternalElkGraphJson.g:1315:5: (lv_properties_1_0= ruleProperty )
                    {
                    // InternalElkGraphJson.g:1315:5: (lv_properties_1_0= ruleProperty )
                    // InternalElkGraphJson.g:1316:6: lv_properties_1_0= ruleProperty
                    {

                    						newCompositeNode(grammarAccess.getElkGraphElementPropertiesAccess().getPropertiesPropertyParserRuleCall_1_0_0());
                    					
                    pushFollow(FOLLOW_4);
                    lv_properties_1_0=ruleProperty();

                    state._fsp--;


                    						if (current==null) {
                    							current = createModelElementForParent(grammarAccess.getElkGraphElementPropertiesRule());
                    						}
                    						add(
                    							current,
                    							"properties",
                    							lv_properties_1_0,
                    							"org.eclipse.elk.graph.json.text.ElkGraphJson.Property");
                    						afterParserOrEnumRuleCall();
                    					

                    }


                    }

                    // InternalElkGraphJson.g:1333:4: (otherlv_2= ',' ( (lv_properties_3_0= ruleProperty ) ) )*
                    loop32:
                    do {
                        int alt32=2;
                        int LA32_0 = input.LA(1);

                        if ( (LA32_0==14) ) {
                            int LA32_1 = input.LA(2);

                            if ( ((LA32_1>=RULE_STRING && LA32_1<=RULE_ID)) ) {
                                alt32=1;
                            }


                        }


                        switch (alt32) {
                    	case 1 :
                    	    // InternalElkGraphJson.g:1334:5: otherlv_2= ',' ( (lv_properties_3_0= ruleProperty ) )
                    	    {
                    	    otherlv_2=(Token)match(input,14,FOLLOW_18); 

                    	    					newLeafNode(otherlv_2, grammarAccess.getElkGraphElementPropertiesAccess().getCommaKeyword_1_1_0());
                    	    				
                    	    // InternalElkGraphJson.g:1338:5: ( (lv_properties_3_0= ruleProperty ) )
                    	    // InternalElkGraphJson.g:1339:6: (lv_properties_3_0= ruleProperty )
                    	    {
                    	    // InternalElkGraphJson.g:1339:6: (lv_properties_3_0= ruleProperty )
                    	    // InternalElkGraphJson.g:1340:7: lv_properties_3_0= ruleProperty
                    	    {

                    	    							newCompositeNode(grammarAccess.getElkGraphElementPropertiesAccess().getPropertiesPropertyParserRuleCall_1_1_1_0());
                    	    						
                    	    pushFollow(FOLLOW_4);
                    	    lv_properties_3_0=ruleProperty();

                    	    state._fsp--;


                    	    							if (current==null) {
                    	    								current = createModelElementForParent(grammarAccess.getElkGraphElementPropertiesRule());
                    	    							}
                    	    							add(
                    	    								current,
                    	    								"properties",
                    	    								lv_properties_3_0,
                    	    								"org.eclipse.elk.graph.json.text.ElkGraphJson.Property");
                    	    							afterParserOrEnumRuleCall();
                    	    						

                    	    }


                    	    }


                    	    }
                    	    break;

                    	default :
                    	    break loop32;
                        }
                    } while (true);


                    }
                    break;

            }

            // InternalElkGraphJson.g:1359:3: (otherlv_4= ',' )?
            int alt34=2;
            int LA34_0 = input.LA(1);

            if ( (LA34_0==14) ) {
                alt34=1;
            }
            switch (alt34) {
                case 1 :
                    // InternalElkGraphJson.g:1360:4: otherlv_4= ','
                    {
                    otherlv_4=(Token)match(input,14,FOLLOW_6); 

                    				newLeafNode(otherlv_4, grammarAccess.getElkGraphElementPropertiesAccess().getCommaKeyword_2());
                    			

                    }
                    break;

            }

            otherlv_5=(Token)match(input,15,FOLLOW_2); 

            			newLeafNode(otherlv_5, grammarAccess.getElkGraphElementPropertiesAccess().getRightCurlyBracketKeyword_3());
            		

            }


            }


            	leaveRule();

        }

            catch (RecognitionException re) {
                recover(input,re);
                appendSkippedTokens();
            }
        finally {
        }
        return current;
    }
    // $ANTLR end "ruleElkGraphElementProperties"


    // $ANTLR start "ruleShapeElement"
    // InternalElkGraphJson.g:1374:1: ruleShapeElement[EObject in_current] returns [EObject current=in_current] : ( ( ruleKeyX otherlv_1= ':' ( (lv_x_2_0= ruleNumber ) ) ) | ( ruleKeyY otherlv_4= ':' ( (lv_y_5_0= ruleNumber ) ) ) | ( ruleKeyWidth otherlv_7= ':' ( (lv_width_8_0= ruleNumber ) ) ) | ( ruleKeyHeight otherlv_10= ':' ( (lv_height_11_0= ruleNumber ) ) ) ) ;
    public final EObject ruleShapeElement(EObject in_current) throws RecognitionException {
        EObject current = in_current;

        Token otherlv_1=null;
        Token otherlv_4=null;
        Token otherlv_7=null;
        Token otherlv_10=null;
        AntlrDatatypeRuleToken lv_x_2_0 = null;

        AntlrDatatypeRuleToken lv_y_5_0 = null;

        AntlrDatatypeRuleToken lv_width_8_0 = null;

        AntlrDatatypeRuleToken lv_height_11_0 = null;



        	enterRule();

        try {
            // InternalElkGraphJson.g:1380:2: ( ( ( ruleKeyX otherlv_1= ':' ( (lv_x_2_0= ruleNumber ) ) ) | ( ruleKeyY otherlv_4= ':' ( (lv_y_5_0= ruleNumber ) ) ) | ( ruleKeyWidth otherlv_7= ':' ( (lv_width_8_0= ruleNumber ) ) ) | ( ruleKeyHeight otherlv_10= ':' ( (lv_height_11_0= ruleNumber ) ) ) ) )
            // InternalElkGraphJson.g:1381:2: ( ( ruleKeyX otherlv_1= ':' ( (lv_x_2_0= ruleNumber ) ) ) | ( ruleKeyY otherlv_4= ':' ( (lv_y_5_0= ruleNumber ) ) ) | ( ruleKeyWidth otherlv_7= ':' ( (lv_width_8_0= ruleNumber ) ) ) | ( ruleKeyHeight otherlv_10= ':' ( (lv_height_11_0= ruleNumber ) ) ) )
            {
            // InternalElkGraphJson.g:1381:2: ( ( ruleKeyX otherlv_1= ':' ( (lv_x_2_0= ruleNumber ) ) ) | ( ruleKeyY otherlv_4= ':' ( (lv_y_5_0= ruleNumber ) ) ) | ( ruleKeyWidth otherlv_7= ':' ( (lv_width_8_0= ruleNumber ) ) ) | ( ruleKeyHeight otherlv_10= ':' ( (lv_height_11_0= ruleNumber ) ) ) )
            int alt35=4;
            switch ( input.LA(1) ) {
            case 43:
            case 44:
            case 45:
                {
                alt35=1;
                }
                break;
            case 46:
            case 47:
            case 48:
                {
                alt35=2;
                }
                break;
            case 49:
            case 50:
            case 51:
                {
                alt35=3;
                }
                break;
            case 52:
            case 53:
            case 54:
                {
                alt35=4;
                }
                break;
            default:
                NoViableAltException nvae =
                    new NoViableAltException("", 35, 0, input);

                throw nvae;
            }

            switch (alt35) {
                case 1 :
                    // InternalElkGraphJson.g:1382:3: ( ruleKeyX otherlv_1= ':' ( (lv_x_2_0= ruleNumber ) ) )
                    {
                    // InternalElkGraphJson.g:1382:3: ( ruleKeyX otherlv_1= ':' ( (lv_x_2_0= ruleNumber ) ) )
                    // InternalElkGraphJson.g:1383:4: ruleKeyX otherlv_1= ':' ( (lv_x_2_0= ruleNumber ) )
                    {

                    				newCompositeNode(grammarAccess.getShapeElementAccess().getKeyXParserRuleCall_0_0());
                    			
                    pushFollow(FOLLOW_7);
                    ruleKeyX();

                    state._fsp--;


                    				afterParserOrEnumRuleCall();
                    			
                    otherlv_1=(Token)match(input,16,FOLLOW_19); 

                    				newLeafNode(otherlv_1, grammarAccess.getShapeElementAccess().getColonKeyword_0_1());
                    			
                    // InternalElkGraphJson.g:1394:4: ( (lv_x_2_0= ruleNumber ) )
                    // InternalElkGraphJson.g:1395:5: (lv_x_2_0= ruleNumber )
                    {
                    // InternalElkGraphJson.g:1395:5: (lv_x_2_0= ruleNumber )
                    // InternalElkGraphJson.g:1396:6: lv_x_2_0= ruleNumber
                    {

                    						newCompositeNode(grammarAccess.getShapeElementAccess().getXNumberParserRuleCall_0_2_0());
                    					
                    pushFollow(FOLLOW_2);
                    lv_x_2_0=ruleNumber();

                    state._fsp--;


                    						if (current==null) {
                    							current = createModelElementForParent(grammarAccess.getShapeElementRule());
                    						}
                    						set(
                    							current,
                    							"x",
                    							lv_x_2_0,
                    							"org.eclipse.elk.graph.json.text.ElkGraphJson.Number");
                    						afterParserOrEnumRuleCall();
                    					

                    }


                    }


                    }


                    }
                    break;
                case 2 :
                    // InternalElkGraphJson.g:1415:3: ( ruleKeyY otherlv_4= ':' ( (lv_y_5_0= ruleNumber ) ) )
                    {
                    // InternalElkGraphJson.g:1415:3: ( ruleKeyY otherlv_4= ':' ( (lv_y_5_0= ruleNumber ) ) )
                    // InternalElkGraphJson.g:1416:4: ruleKeyY otherlv_4= ':' ( (lv_y_5_0= ruleNumber ) )
                    {

                    				newCompositeNode(grammarAccess.getShapeElementAccess().getKeyYParserRuleCall_1_0());
                    			
                    pushFollow(FOLLOW_7);
                    ruleKeyY();

                    state._fsp--;


                    				afterParserOrEnumRuleCall();
                    			
                    otherlv_4=(Token)match(input,16,FOLLOW_19); 

                    				newLeafNode(otherlv_4, grammarAccess.getShapeElementAccess().getColonKeyword_1_1());
                    			
                    // InternalElkGraphJson.g:1427:4: ( (lv_y_5_0= ruleNumber ) )
                    // InternalElkGraphJson.g:1428:5: (lv_y_5_0= ruleNumber )
                    {
                    // InternalElkGraphJson.g:1428:5: (lv_y_5_0= ruleNumber )
                    // InternalElkGraphJson.g:1429:6: lv_y_5_0= ruleNumber
                    {

                    						newCompositeNode(grammarAccess.getShapeElementAccess().getYNumberParserRuleCall_1_2_0());
                    					
                    pushFollow(FOLLOW_2);
                    lv_y_5_0=ruleNumber();

                    state._fsp--;


                    						if (current==null) {
                    							current = createModelElementForParent(grammarAccess.getShapeElementRule());
                    						}
                    						set(
                    							current,
                    							"y",
                    							lv_y_5_0,
                    							"org.eclipse.elk.graph.json.text.ElkGraphJson.Number");
                    						afterParserOrEnumRuleCall();
                    					

                    }


                    }


                    }


                    }
                    break;
                case 3 :
                    // InternalElkGraphJson.g:1448:3: ( ruleKeyWidth otherlv_7= ':' ( (lv_width_8_0= ruleNumber ) ) )
                    {
                    // InternalElkGraphJson.g:1448:3: ( ruleKeyWidth otherlv_7= ':' ( (lv_width_8_0= ruleNumber ) ) )
                    // InternalElkGraphJson.g:1449:4: ruleKeyWidth otherlv_7= ':' ( (lv_width_8_0= ruleNumber ) )
                    {

                    				newCompositeNode(grammarAccess.getShapeElementAccess().getKeyWidthParserRuleCall_2_0());
                    			
                    pushFollow(FOLLOW_7);
                    ruleKeyWidth();

                    state._fsp--;


                    				afterParserOrEnumRuleCall();
                    			
                    otherlv_7=(Token)match(input,16,FOLLOW_19); 

                    				newLeafNode(otherlv_7, grammarAccess.getShapeElementAccess().getColonKeyword_2_1());
                    			
                    // InternalElkGraphJson.g:1460:4: ( (lv_width_8_0= ruleNumber ) )
                    // InternalElkGraphJson.g:1461:5: (lv_width_8_0= ruleNumber )
                    {
                    // InternalElkGraphJson.g:1461:5: (lv_width_8_0= ruleNumber )
                    // InternalElkGraphJson.g:1462:6: lv_width_8_0= ruleNumber
                    {

                    						newCompositeNode(grammarAccess.getShapeElementAccess().getWidthNumberParserRuleCall_2_2_0());
                    					
                    pushFollow(FOLLOW_2);
                    lv_width_8_0=ruleNumber();

                    state._fsp--;


                    						if (current==null) {
                    							current = createModelElementForParent(grammarAccess.getShapeElementRule());
                    						}
                    						set(
                    							current,
                    							"width",
                    							lv_width_8_0,
                    							"org.eclipse.elk.graph.json.text.ElkGraphJson.Number");
                    						afterParserOrEnumRuleCall();
                    					

                    }


                    }


                    }


                    }
                    break;
                case 4 :
                    // InternalElkGraphJson.g:1481:3: ( ruleKeyHeight otherlv_10= ':' ( (lv_height_11_0= ruleNumber ) ) )
                    {
                    // InternalElkGraphJson.g:1481:3: ( ruleKeyHeight otherlv_10= ':' ( (lv_height_11_0= ruleNumber ) ) )
                    // InternalElkGraphJson.g:1482:4: ruleKeyHeight otherlv_10= ':' ( (lv_height_11_0= ruleNumber ) )
                    {

                    				newCompositeNode(grammarAccess.getShapeElementAccess().getKeyHeightParserRuleCall_3_0());
                    			
                    pushFollow(FOLLOW_7);
                    ruleKeyHeight();

                    state._fsp--;


                    				afterParserOrEnumRuleCall();
                    			
                    otherlv_10=(Token)match(input,16,FOLLOW_19); 

                    				newLeafNode(otherlv_10, grammarAccess.getShapeElementAccess().getColonKeyword_3_1());
                    			
                    // InternalElkGraphJson.g:1493:4: ( (lv_height_11_0= ruleNumber ) )
                    // InternalElkGraphJson.g:1494:5: (lv_height_11_0= ruleNumber )
                    {
                    // InternalElkGraphJson.g:1494:5: (lv_height_11_0= ruleNumber )
                    // InternalElkGraphJson.g:1495:6: lv_height_11_0= ruleNumber
                    {

                    						newCompositeNode(grammarAccess.getShapeElementAccess().getHeightNumberParserRuleCall_3_2_0());
                    					
                    pushFollow(FOLLOW_2);
                    lv_height_11_0=ruleNumber();

                    state._fsp--;


                    						if (current==null) {
                    							current = createModelElementForParent(grammarAccess.getShapeElementRule());
                    						}
                    						set(
                    							current,
                    							"height",
                    							lv_height_11_0,
                    							"org.eclipse.elk.graph.json.text.ElkGraphJson.Number");
                    						afterParserOrEnumRuleCall();
                    					

                    }


                    }


                    }


                    }
                    break;

            }


            }


            	leaveRule();

        }

            catch (RecognitionException re) {
                recover(input,re);
                appendSkippedTokens();
            }
        finally {
        }
        return current;
    }
    // $ANTLR end "ruleShapeElement"


    // $ANTLR start "entryRuleProperty"
    // InternalElkGraphJson.g:1517:1: entryRuleProperty returns [EObject current=null] : iv_ruleProperty= ruleProperty EOF ;
    public final EObject entryRuleProperty() throws RecognitionException {
        EObject current = null;

        EObject iv_ruleProperty = null;


        try {
            // InternalElkGraphJson.g:1517:49: (iv_ruleProperty= ruleProperty EOF )
            // InternalElkGraphJson.g:1518:2: iv_ruleProperty= ruleProperty EOF
            {
             newCompositeNode(grammarAccess.getPropertyRule()); 
            pushFollow(FOLLOW_1);
            iv_ruleProperty=ruleProperty();

            state._fsp--;

             current =iv_ruleProperty; 
            match(input,EOF,FOLLOW_2); 

            }

        }

            catch (RecognitionException re) {
                recover(input,re);
                appendSkippedTokens();
            }
        finally {
        }
        return current;
    }
    // $ANTLR end "entryRuleProperty"


    // $ANTLR start "ruleProperty"
    // InternalElkGraphJson.g:1524:1: ruleProperty returns [EObject current=null] : ( ( (lv_key_0_0= rulePropertyKey ) ) otherlv_1= ':' ( ( (lv_value_2_0= ruleStringValue ) ) | ( (lv_value_3_0= ruleNumberValue ) ) | ( (lv_value_4_0= ruleBooleanValue ) ) ) ) ;
    public final EObject ruleProperty() throws RecognitionException {
        EObject current = null;

        Token otherlv_1=null;
        AntlrDatatypeRuleToken lv_key_0_0 = null;

        AntlrDatatypeRuleToken lv_value_2_0 = null;

        AntlrDatatypeRuleToken lv_value_3_0 = null;

        AntlrDatatypeRuleToken lv_value_4_0 = null;



        	enterRule();

        try {
            // InternalElkGraphJson.g:1530:2: ( ( ( (lv_key_0_0= rulePropertyKey ) ) otherlv_1= ':' ( ( (lv_value_2_0= ruleStringValue ) ) | ( (lv_value_3_0= ruleNumberValue ) ) | ( (lv_value_4_0= ruleBooleanValue ) ) ) ) )
            // InternalElkGraphJson.g:1531:2: ( ( (lv_key_0_0= rulePropertyKey ) ) otherlv_1= ':' ( ( (lv_value_2_0= ruleStringValue ) ) | ( (lv_value_3_0= ruleNumberValue ) ) | ( (lv_value_4_0= ruleBooleanValue ) ) ) )
            {
            // InternalElkGraphJson.g:1531:2: ( ( (lv_key_0_0= rulePropertyKey ) ) otherlv_1= ':' ( ( (lv_value_2_0= ruleStringValue ) ) | ( (lv_value_3_0= ruleNumberValue ) ) | ( (lv_value_4_0= ruleBooleanValue ) ) ) )
            // InternalElkGraphJson.g:1532:3: ( (lv_key_0_0= rulePropertyKey ) ) otherlv_1= ':' ( ( (lv_value_2_0= ruleStringValue ) ) | ( (lv_value_3_0= ruleNumberValue ) ) | ( (lv_value_4_0= ruleBooleanValue ) ) )
            {
            // InternalElkGraphJson.g:1532:3: ( (lv_key_0_0= rulePropertyKey ) )
            // InternalElkGraphJson.g:1533:4: (lv_key_0_0= rulePropertyKey )
            {
            // InternalElkGraphJson.g:1533:4: (lv_key_0_0= rulePropertyKey )
            // InternalElkGraphJson.g:1534:5: lv_key_0_0= rulePropertyKey
            {

            					newCompositeNode(grammarAccess.getPropertyAccess().getKeyPropertyKeyParserRuleCall_0_0());
            				
            pushFollow(FOLLOW_7);
            lv_key_0_0=rulePropertyKey();

            state._fsp--;


            					if (current==null) {
            						current = createModelElementForParent(grammarAccess.getPropertyRule());
            					}
            					set(
            						current,
            						"key",
            						lv_key_0_0,
            						"org.eclipse.elk.graph.json.text.ElkGraphJson.PropertyKey");
            					afterParserOrEnumRuleCall();
            				

            }


            }

            otherlv_1=(Token)match(input,16,FOLLOW_20); 

            			newLeafNode(otherlv_1, grammarAccess.getPropertyAccess().getColonKeyword_1());
            		
            // InternalElkGraphJson.g:1555:3: ( ( (lv_value_2_0= ruleStringValue ) ) | ( (lv_value_3_0= ruleNumberValue ) ) | ( (lv_value_4_0= ruleBooleanValue ) ) )
            int alt36=3;
            switch ( input.LA(1) ) {
            case RULE_STRING:
                {
                alt36=1;
                }
                break;
            case RULE_SIGNED_INT:
            case RULE_FLOAT:
                {
                alt36=2;
                }
                break;
            case 19:
            case 20:
                {
                alt36=3;
                }
                break;
            default:
                NoViableAltException nvae =
                    new NoViableAltException("", 36, 0, input);

                throw nvae;
            }

            switch (alt36) {
                case 1 :
                    // InternalElkGraphJson.g:1556:4: ( (lv_value_2_0= ruleStringValue ) )
                    {
                    // InternalElkGraphJson.g:1556:4: ( (lv_value_2_0= ruleStringValue ) )
                    // InternalElkGraphJson.g:1557:5: (lv_value_2_0= ruleStringValue )
                    {
                    // InternalElkGraphJson.g:1557:5: (lv_value_2_0= ruleStringValue )
                    // InternalElkGraphJson.g:1558:6: lv_value_2_0= ruleStringValue
                    {

                    						newCompositeNode(grammarAccess.getPropertyAccess().getValueStringValueParserRuleCall_2_0_0());
                    					
                    pushFollow(FOLLOW_2);
                    lv_value_2_0=ruleStringValue();

                    state._fsp--;


                    						if (current==null) {
                    							current = createModelElementForParent(grammarAccess.getPropertyRule());
                    						}
                    						set(
                    							current,
                    							"value",
                    							lv_value_2_0,
                    							"org.eclipse.elk.graph.json.text.ElkGraphJson.StringValue");
                    						afterParserOrEnumRuleCall();
                    					

                    }


                    }


                    }
                    break;
                case 2 :
                    // InternalElkGraphJson.g:1576:4: ( (lv_value_3_0= ruleNumberValue ) )
                    {
                    // InternalElkGraphJson.g:1576:4: ( (lv_value_3_0= ruleNumberValue ) )
                    // InternalElkGraphJson.g:1577:5: (lv_value_3_0= ruleNumberValue )
                    {
                    // InternalElkGraphJson.g:1577:5: (lv_value_3_0= ruleNumberValue )
                    // InternalElkGraphJson.g:1578:6: lv_value_3_0= ruleNumberValue
                    {

                    						newCompositeNode(grammarAccess.getPropertyAccess().getValueNumberValueParserRuleCall_2_1_0());
                    					
                    pushFollow(FOLLOW_2);
                    lv_value_3_0=ruleNumberValue();

                    state._fsp--;


                    						if (current==null) {
                    							current = createModelElementForParent(grammarAccess.getPropertyRule());
                    						}
                    						set(
                    							current,
                    							"value",
                    							lv_value_3_0,
                    							"org.eclipse.elk.graph.json.text.ElkGraphJson.NumberValue");
                    						afterParserOrEnumRuleCall();
                    					

                    }


                    }


                    }
                    break;
                case 3 :
                    // InternalElkGraphJson.g:1596:4: ( (lv_value_4_0= ruleBooleanValue ) )
                    {
                    // InternalElkGraphJson.g:1596:4: ( (lv_value_4_0= ruleBooleanValue ) )
                    // InternalElkGraphJson.g:1597:5: (lv_value_4_0= ruleBooleanValue )
                    {
                    // InternalElkGraphJson.g:1597:5: (lv_value_4_0= ruleBooleanValue )
                    // InternalElkGraphJson.g:1598:6: lv_value_4_0= ruleBooleanValue
                    {

                    						newCompositeNode(grammarAccess.getPropertyAccess().getValueBooleanValueParserRuleCall_2_2_0());
                    					
                    pushFollow(FOLLOW_2);
                    lv_value_4_0=ruleBooleanValue();

                    state._fsp--;


                    						if (current==null) {
                    							current = createModelElementForParent(grammarAccess.getPropertyRule());
                    						}
                    						set(
                    							current,
                    							"value",
                    							lv_value_4_0,
                    							"org.eclipse.elk.graph.json.text.ElkGraphJson.BooleanValue");
                    						afterParserOrEnumRuleCall();
                    					

                    }


                    }


                    }
                    break;

            }


            }


            }


            	leaveRule();

        }

            catch (RecognitionException re) {
                recover(input,re);
                appendSkippedTokens();
            }
        finally {
        }
        return current;
    }
    // $ANTLR end "ruleProperty"


    // $ANTLR start "entryRulePropertyKey"
    // InternalElkGraphJson.g:1620:1: entryRulePropertyKey returns [String current=null] : iv_rulePropertyKey= rulePropertyKey EOF ;
    public final String entryRulePropertyKey() throws RecognitionException {
        String current = null;

        AntlrDatatypeRuleToken iv_rulePropertyKey = null;



        	HiddenTokens myHiddenTokenState = ((XtextTokenStream)input).setHiddenTokens();

        try {
            // InternalElkGraphJson.g:1622:2: (iv_rulePropertyKey= rulePropertyKey EOF )
            // InternalElkGraphJson.g:1623:2: iv_rulePropertyKey= rulePropertyKey EOF
            {
             newCompositeNode(grammarAccess.getPropertyKeyRule()); 
            pushFollow(FOLLOW_1);
            iv_rulePropertyKey=rulePropertyKey();

            state._fsp--;

             current =iv_rulePropertyKey.getText(); 
            match(input,EOF,FOLLOW_2); 

            }

        }

            catch (RecognitionException re) {
                recover(input,re);
                appendSkippedTokens();
            }
        finally {

            	myHiddenTokenState.restore();

        }
        return current;
    }
    // $ANTLR end "entryRulePropertyKey"


    // $ANTLR start "rulePropertyKey"
    // InternalElkGraphJson.g:1632:1: rulePropertyKey returns [AntlrDatatypeRuleToken current=new AntlrDatatypeRuleToken()] : (this_STRING_0= RULE_STRING | this_ID_1= RULE_ID ) ;
    public final AntlrDatatypeRuleToken rulePropertyKey() throws RecognitionException {
        AntlrDatatypeRuleToken current = new AntlrDatatypeRuleToken();

        Token this_STRING_0=null;
        Token this_ID_1=null;


        	enterRule();
        	HiddenTokens myHiddenTokenState = ((XtextTokenStream)input).setHiddenTokens();

        try {
            // InternalElkGraphJson.g:1639:2: ( (this_STRING_0= RULE_STRING | this_ID_1= RULE_ID ) )
            // InternalElkGraphJson.g:1640:2: (this_STRING_0= RULE_STRING | this_ID_1= RULE_ID )
            {
            // InternalElkGraphJson.g:1640:2: (this_STRING_0= RULE_STRING | this_ID_1= RULE_ID )
            int alt37=2;
            int LA37_0 = input.LA(1);

            if ( (LA37_0==RULE_STRING) ) {
                alt37=1;
            }
            else if ( (LA37_0==RULE_ID) ) {
                alt37=2;
            }
            else {
                NoViableAltException nvae =
                    new NoViableAltException("", 37, 0, input);

                throw nvae;
            }
            switch (alt37) {
                case 1 :
                    // InternalElkGraphJson.g:1641:3: this_STRING_0= RULE_STRING
                    {
                    this_STRING_0=(Token)match(input,RULE_STRING,FOLLOW_2); 

                    			current.merge(this_STRING_0);
                    		

                    			newLeafNode(this_STRING_0, grammarAccess.getPropertyKeyAccess().getSTRINGTerminalRuleCall_0());
                    		

                    }
                    break;
                case 2 :
                    // InternalElkGraphJson.g:1649:3: this_ID_1= RULE_ID
                    {
                    this_ID_1=(Token)match(input,RULE_ID,FOLLOW_2); 

                    			current.merge(this_ID_1);
                    		

                    			newLeafNode(this_ID_1, grammarAccess.getPropertyKeyAccess().getIDTerminalRuleCall_1());
                    		

                    }
                    break;

            }


            }


            	leaveRule();

        }

            catch (RecognitionException re) {
                recover(input,re);
                appendSkippedTokens();
            }
        finally {

            	myHiddenTokenState.restore();

        }
        return current;
    }
    // $ANTLR end "rulePropertyKey"


    // $ANTLR start "entryRuleStringValue"
    // InternalElkGraphJson.g:1663:1: entryRuleStringValue returns [String current=null] : iv_ruleStringValue= ruleStringValue EOF ;
    public final String entryRuleStringValue() throws RecognitionException {
        String current = null;

        AntlrDatatypeRuleToken iv_ruleStringValue = null;


        try {
            // InternalElkGraphJson.g:1663:51: (iv_ruleStringValue= ruleStringValue EOF )
            // InternalElkGraphJson.g:1664:2: iv_ruleStringValue= ruleStringValue EOF
            {
             newCompositeNode(grammarAccess.getStringValueRule()); 
            pushFollow(FOLLOW_1);
            iv_ruleStringValue=ruleStringValue();

            state._fsp--;

             current =iv_ruleStringValue.getText(); 
            match(input,EOF,FOLLOW_2); 

            }

        }

            catch (RecognitionException re) {
                recover(input,re);
                appendSkippedTokens();
            }
        finally {
        }
        return current;
    }
    // $ANTLR end "entryRuleStringValue"


    // $ANTLR start "ruleStringValue"
    // InternalElkGraphJson.g:1670:1: ruleStringValue returns [AntlrDatatypeRuleToken current=new AntlrDatatypeRuleToken()] : this_STRING_0= RULE_STRING ;
    public final AntlrDatatypeRuleToken ruleStringValue() throws RecognitionException {
        AntlrDatatypeRuleToken current = new AntlrDatatypeRuleToken();

        Token this_STRING_0=null;


        	enterRule();

        try {
            // InternalElkGraphJson.g:1676:2: (this_STRING_0= RULE_STRING )
            // InternalElkGraphJson.g:1677:2: this_STRING_0= RULE_STRING
            {
            this_STRING_0=(Token)match(input,RULE_STRING,FOLLOW_2); 

            		current.merge(this_STRING_0);
            	

            		newLeafNode(this_STRING_0, grammarAccess.getStringValueAccess().getSTRINGTerminalRuleCall());
            	

            }


            	leaveRule();

        }

            catch (RecognitionException re) {
                recover(input,re);
                appendSkippedTokens();
            }
        finally {
        }
        return current;
    }
    // $ANTLR end "ruleStringValue"


    // $ANTLR start "entryRuleNumberValue"
    // InternalElkGraphJson.g:1687:1: entryRuleNumberValue returns [String current=null] : iv_ruleNumberValue= ruleNumberValue EOF ;
    public final String entryRuleNumberValue() throws RecognitionException {
        String current = null;

        AntlrDatatypeRuleToken iv_ruleNumberValue = null;


        try {
            // InternalElkGraphJson.g:1687:51: (iv_ruleNumberValue= ruleNumberValue EOF )
            // InternalElkGraphJson.g:1688:2: iv_ruleNumberValue= ruleNumberValue EOF
            {
             newCompositeNode(grammarAccess.getNumberValueRule()); 
            pushFollow(FOLLOW_1);
            iv_ruleNumberValue=ruleNumberValue();

            state._fsp--;

             current =iv_ruleNumberValue.getText(); 
            match(input,EOF,FOLLOW_2); 

            }

        }

            catch (RecognitionException re) {
                recover(input,re);
                appendSkippedTokens();
            }
        finally {
        }
        return current;
    }
    // $ANTLR end "entryRuleNumberValue"


    // $ANTLR start "ruleNumberValue"
    // InternalElkGraphJson.g:1694:1: ruleNumberValue returns [AntlrDatatypeRuleToken current=new AntlrDatatypeRuleToken()] : (this_SIGNED_INT_0= RULE_SIGNED_INT | this_FLOAT_1= RULE_FLOAT ) ;
    public final AntlrDatatypeRuleToken ruleNumberValue() throws RecognitionException {
        AntlrDatatypeRuleToken current = new AntlrDatatypeRuleToken();

        Token this_SIGNED_INT_0=null;
        Token this_FLOAT_1=null;


        	enterRule();

        try {
            // InternalElkGraphJson.g:1700:2: ( (this_SIGNED_INT_0= RULE_SIGNED_INT | this_FLOAT_1= RULE_FLOAT ) )
            // InternalElkGraphJson.g:1701:2: (this_SIGNED_INT_0= RULE_SIGNED_INT | this_FLOAT_1= RULE_FLOAT )
            {
            // InternalElkGraphJson.g:1701:2: (this_SIGNED_INT_0= RULE_SIGNED_INT | this_FLOAT_1= RULE_FLOAT )
            int alt38=2;
            int LA38_0 = input.LA(1);

            if ( (LA38_0==RULE_SIGNED_INT) ) {
                alt38=1;
            }
            else if ( (LA38_0==RULE_FLOAT) ) {
                alt38=2;
            }
            else {
                NoViableAltException nvae =
                    new NoViableAltException("", 38, 0, input);

                throw nvae;
            }
            switch (alt38) {
                case 1 :
                    // InternalElkGraphJson.g:1702:3: this_SIGNED_INT_0= RULE_SIGNED_INT
                    {
                    this_SIGNED_INT_0=(Token)match(input,RULE_SIGNED_INT,FOLLOW_2); 

                    			current.merge(this_SIGNED_INT_0);
                    		

                    			newLeafNode(this_SIGNED_INT_0, grammarAccess.getNumberValueAccess().getSIGNED_INTTerminalRuleCall_0());
                    		

                    }
                    break;
                case 2 :
                    // InternalElkGraphJson.g:1710:3: this_FLOAT_1= RULE_FLOAT
                    {
                    this_FLOAT_1=(Token)match(input,RULE_FLOAT,FOLLOW_2); 

                    			current.merge(this_FLOAT_1);
                    		

                    			newLeafNode(this_FLOAT_1, grammarAccess.getNumberValueAccess().getFLOATTerminalRuleCall_1());
                    		

                    }
                    break;

            }


            }


            	leaveRule();

        }

            catch (RecognitionException re) {
                recover(input,re);
                appendSkippedTokens();
            }
        finally {
        }
        return current;
    }
    // $ANTLR end "ruleNumberValue"


    // $ANTLR start "entryRuleBooleanValue"
    // InternalElkGraphJson.g:1721:1: entryRuleBooleanValue returns [String current=null] : iv_ruleBooleanValue= ruleBooleanValue EOF ;
    public final String entryRuleBooleanValue() throws RecognitionException {
        String current = null;

        AntlrDatatypeRuleToken iv_ruleBooleanValue = null;


        try {
            // InternalElkGraphJson.g:1721:52: (iv_ruleBooleanValue= ruleBooleanValue EOF )
            // InternalElkGraphJson.g:1722:2: iv_ruleBooleanValue= ruleBooleanValue EOF
            {
             newCompositeNode(grammarAccess.getBooleanValueRule()); 
            pushFollow(FOLLOW_1);
            iv_ruleBooleanValue=ruleBooleanValue();

            state._fsp--;

             current =iv_ruleBooleanValue.getText(); 
            match(input,EOF,FOLLOW_2); 

            }

        }

            catch (RecognitionException re) {
                recover(input,re);
                appendSkippedTokens();
            }
        finally {
        }
        return current;
    }
    // $ANTLR end "entryRuleBooleanValue"


    // $ANTLR start "ruleBooleanValue"
    // InternalElkGraphJson.g:1728:1: ruleBooleanValue returns [AntlrDatatypeRuleToken current=new AntlrDatatypeRuleToken()] : (kw= 'true' | kw= 'false' ) ;
    public final AntlrDatatypeRuleToken ruleBooleanValue() throws RecognitionException {
        AntlrDatatypeRuleToken current = new AntlrDatatypeRuleToken();

        Token kw=null;


        	enterRule();

        try {
            // InternalElkGraphJson.g:1734:2: ( (kw= 'true' | kw= 'false' ) )
            // InternalElkGraphJson.g:1735:2: (kw= 'true' | kw= 'false' )
            {
            // InternalElkGraphJson.g:1735:2: (kw= 'true' | kw= 'false' )
            int alt39=2;
            int LA39_0 = input.LA(1);

            if ( (LA39_0==19) ) {
                alt39=1;
            }
            else if ( (LA39_0==20) ) {
                alt39=2;
            }
            else {
                NoViableAltException nvae =
                    new NoViableAltException("", 39, 0, input);

                throw nvae;
            }
            switch (alt39) {
                case 1 :
                    // InternalElkGraphJson.g:1736:3: kw= 'true'
                    {
                    kw=(Token)match(input,19,FOLLOW_2); 

                    			current.merge(kw);
                    			newLeafNode(kw, grammarAccess.getBooleanValueAccess().getTrueKeyword_0());
                    		

                    }
                    break;
                case 2 :
                    // InternalElkGraphJson.g:1742:3: kw= 'false'
                    {
                    kw=(Token)match(input,20,FOLLOW_2); 

                    			current.merge(kw);
                    			newLeafNode(kw, grammarAccess.getBooleanValueAccess().getFalseKeyword_1());
                    		

                    }
                    break;

            }


            }


            	leaveRule();

        }

            catch (RecognitionException re) {
                recover(input,re);
                appendSkippedTokens();
            }
        finally {
        }
        return current;
    }
    // $ANTLR end "ruleBooleanValue"


    // $ANTLR start "entryRuleNumber"
    // InternalElkGraphJson.g:1751:1: entryRuleNumber returns [String current=null] : iv_ruleNumber= ruleNumber EOF ;
    public final String entryRuleNumber() throws RecognitionException {
        String current = null;

        AntlrDatatypeRuleToken iv_ruleNumber = null;


        try {
            // InternalElkGraphJson.g:1751:46: (iv_ruleNumber= ruleNumber EOF )
            // InternalElkGraphJson.g:1752:2: iv_ruleNumber= ruleNumber EOF
            {
             newCompositeNode(grammarAccess.getNumberRule()); 
            pushFollow(FOLLOW_1);
            iv_ruleNumber=ruleNumber();

            state._fsp--;

             current =iv_ruleNumber.getText(); 
            match(input,EOF,FOLLOW_2); 

            }

        }

            catch (RecognitionException re) {
                recover(input,re);
                appendSkippedTokens();
            }
        finally {
        }
        return current;
    }
    // $ANTLR end "entryRuleNumber"


    // $ANTLR start "ruleNumber"
    // InternalElkGraphJson.g:1758:1: ruleNumber returns [AntlrDatatypeRuleToken current=new AntlrDatatypeRuleToken()] : (this_SIGNED_INT_0= RULE_SIGNED_INT | this_FLOAT_1= RULE_FLOAT ) ;
    public final AntlrDatatypeRuleToken ruleNumber() throws RecognitionException {
        AntlrDatatypeRuleToken current = new AntlrDatatypeRuleToken();

        Token this_SIGNED_INT_0=null;
        Token this_FLOAT_1=null;


        	enterRule();

        try {
            // InternalElkGraphJson.g:1764:2: ( (this_SIGNED_INT_0= RULE_SIGNED_INT | this_FLOAT_1= RULE_FLOAT ) )
            // InternalElkGraphJson.g:1765:2: (this_SIGNED_INT_0= RULE_SIGNED_INT | this_FLOAT_1= RULE_FLOAT )
            {
            // InternalElkGraphJson.g:1765:2: (this_SIGNED_INT_0= RULE_SIGNED_INT | this_FLOAT_1= RULE_FLOAT )
            int alt40=2;
            int LA40_0 = input.LA(1);

            if ( (LA40_0==RULE_SIGNED_INT) ) {
                alt40=1;
            }
            else if ( (LA40_0==RULE_FLOAT) ) {
                alt40=2;
            }
            else {
                NoViableAltException nvae =
                    new NoViableAltException("", 40, 0, input);

                throw nvae;
            }
            switch (alt40) {
                case 1 :
                    // InternalElkGraphJson.g:1766:3: this_SIGNED_INT_0= RULE_SIGNED_INT
                    {
                    this_SIGNED_INT_0=(Token)match(input,RULE_SIGNED_INT,FOLLOW_2); 

                    			current.merge(this_SIGNED_INT_0);
                    		

                    			newLeafNode(this_SIGNED_INT_0, grammarAccess.getNumberAccess().getSIGNED_INTTerminalRuleCall_0());
                    		

                    }
                    break;
                case 2 :
                    // InternalElkGraphJson.g:1774:3: this_FLOAT_1= RULE_FLOAT
                    {
                    this_FLOAT_1=(Token)match(input,RULE_FLOAT,FOLLOW_2); 

                    			current.merge(this_FLOAT_1);
                    		

                    			newLeafNode(this_FLOAT_1, grammarAccess.getNumberAccess().getFLOATTerminalRuleCall_1());
                    		

                    }
                    break;

            }


            }


            	leaveRule();

        }

            catch (RecognitionException re) {
                recover(input,re);
                appendSkippedTokens();
            }
        finally {
        }
        return current;
    }
    // $ANTLR end "ruleNumber"


    // $ANTLR start "entryRuleJsonObject"
    // InternalElkGraphJson.g:1785:1: entryRuleJsonObject returns [String current=null] : iv_ruleJsonObject= ruleJsonObject EOF ;
    public final String entryRuleJsonObject() throws RecognitionException {
        String current = null;

        AntlrDatatypeRuleToken iv_ruleJsonObject = null;


        try {
            // InternalElkGraphJson.g:1785:50: (iv_ruleJsonObject= ruleJsonObject EOF )
            // InternalElkGraphJson.g:1786:2: iv_ruleJsonObject= ruleJsonObject EOF
            {
             newCompositeNode(grammarAccess.getJsonObjectRule()); 
            pushFollow(FOLLOW_1);
            iv_ruleJsonObject=ruleJsonObject();

            state._fsp--;

             current =iv_ruleJsonObject.getText(); 
            match(input,EOF,FOLLOW_2); 

            }

        }

            catch (RecognitionException re) {
                recover(input,re);
                appendSkippedTokens();
            }
        finally {
        }
        return current;
    }
    // $ANTLR end "entryRuleJsonObject"


    // $ANTLR start "ruleJsonObject"
    // InternalElkGraphJson.g:1792:1: ruleJsonObject returns [AntlrDatatypeRuleToken current=new AntlrDatatypeRuleToken()] : (kw= '{' (this_JsonMember_1= ruleJsonMember (kw= ',' this_JsonMember_3= ruleJsonMember )* )? (kw= ',' )? kw= '}' ) ;
    public final AntlrDatatypeRuleToken ruleJsonObject() throws RecognitionException {
        AntlrDatatypeRuleToken current = new AntlrDatatypeRuleToken();

        Token kw=null;
        AntlrDatatypeRuleToken this_JsonMember_1 = null;

        AntlrDatatypeRuleToken this_JsonMember_3 = null;



        	enterRule();

        try {
            // InternalElkGraphJson.g:1798:2: ( (kw= '{' (this_JsonMember_1= ruleJsonMember (kw= ',' this_JsonMember_3= ruleJsonMember )* )? (kw= ',' )? kw= '}' ) )
            // InternalElkGraphJson.g:1799:2: (kw= '{' (this_JsonMember_1= ruleJsonMember (kw= ',' this_JsonMember_3= ruleJsonMember )* )? (kw= ',' )? kw= '}' )
            {
            // InternalElkGraphJson.g:1799:2: (kw= '{' (this_JsonMember_1= ruleJsonMember (kw= ',' this_JsonMember_3= ruleJsonMember )* )? (kw= ',' )? kw= '}' )
            // InternalElkGraphJson.g:1800:3: kw= '{' (this_JsonMember_1= ruleJsonMember (kw= ',' this_JsonMember_3= ruleJsonMember )* )? (kw= ',' )? kw= '}'
            {
            kw=(Token)match(input,13,FOLLOW_3); 

            			current.merge(kw);
            			newLeafNode(kw, grammarAccess.getJsonObjectAccess().getLeftCurlyBracketKeyword_0());
            		
            // InternalElkGraphJson.g:1805:3: (this_JsonMember_1= ruleJsonMember (kw= ',' this_JsonMember_3= ruleJsonMember )* )?
            int alt42=2;
            int LA42_0 = input.LA(1);

            if ( ((LA42_0>=RULE_STRING && LA42_0<=RULE_ID)) ) {
                alt42=1;
            }
            switch (alt42) {
                case 1 :
                    // InternalElkGraphJson.g:1806:4: this_JsonMember_1= ruleJsonMember (kw= ',' this_JsonMember_3= ruleJsonMember )*
                    {

                    				newCompositeNode(grammarAccess.getJsonObjectAccess().getJsonMemberParserRuleCall_1_0());
                    			
                    pushFollow(FOLLOW_4);
                    this_JsonMember_1=ruleJsonMember();

                    state._fsp--;


                    				current.merge(this_JsonMember_1);
                    			

                    				afterParserOrEnumRuleCall();
                    			
                    // InternalElkGraphJson.g:1816:4: (kw= ',' this_JsonMember_3= ruleJsonMember )*
                    loop41:
                    do {
                        int alt41=2;
                        int LA41_0 = input.LA(1);

                        if ( (LA41_0==14) ) {
                            int LA41_1 = input.LA(2);

                            if ( ((LA41_1>=RULE_STRING && LA41_1<=RULE_ID)) ) {
                                alt41=1;
                            }


                        }


                        switch (alt41) {
                    	case 1 :
                    	    // InternalElkGraphJson.g:1817:5: kw= ',' this_JsonMember_3= ruleJsonMember
                    	    {
                    	    kw=(Token)match(input,14,FOLLOW_5); 

                    	    					current.merge(kw);
                    	    					newLeafNode(kw, grammarAccess.getJsonObjectAccess().getCommaKeyword_1_1_0());
                    	    				

                    	    					newCompositeNode(grammarAccess.getJsonObjectAccess().getJsonMemberParserRuleCall_1_1_1());
                    	    				
                    	    pushFollow(FOLLOW_4);
                    	    this_JsonMember_3=ruleJsonMember();

                    	    state._fsp--;


                    	    					current.merge(this_JsonMember_3);
                    	    				

                    	    					afterParserOrEnumRuleCall();
                    	    				

                    	    }
                    	    break;

                    	default :
                    	    break loop41;
                        }
                    } while (true);


                    }
                    break;

            }

            // InternalElkGraphJson.g:1834:3: (kw= ',' )?
            int alt43=2;
            int LA43_0 = input.LA(1);

            if ( (LA43_0==14) ) {
                alt43=1;
            }
            switch (alt43) {
                case 1 :
                    // InternalElkGraphJson.g:1835:4: kw= ','
                    {
                    kw=(Token)match(input,14,FOLLOW_6); 

                    				current.merge(kw);
                    				newLeafNode(kw, grammarAccess.getJsonObjectAccess().getCommaKeyword_2());
                    			

                    }
                    break;

            }

            kw=(Token)match(input,15,FOLLOW_2); 

            			current.merge(kw);
            			newLeafNode(kw, grammarAccess.getJsonObjectAccess().getRightCurlyBracketKeyword_3());
            		

            }


            }


            	leaveRule();

        }

            catch (RecognitionException re) {
                recover(input,re);
                appendSkippedTokens();
            }
        finally {
        }
        return current;
    }
    // $ANTLR end "ruleJsonObject"


    // $ANTLR start "entryRuleJsonArray"
    // InternalElkGraphJson.g:1850:1: entryRuleJsonArray returns [String current=null] : iv_ruleJsonArray= ruleJsonArray EOF ;
    public final String entryRuleJsonArray() throws RecognitionException {
        String current = null;

        AntlrDatatypeRuleToken iv_ruleJsonArray = null;


        try {
            // InternalElkGraphJson.g:1850:49: (iv_ruleJsonArray= ruleJsonArray EOF )
            // InternalElkGraphJson.g:1851:2: iv_ruleJsonArray= ruleJsonArray EOF
            {
             newCompositeNode(grammarAccess.getJsonArrayRule()); 
            pushFollow(FOLLOW_1);
            iv_ruleJsonArray=ruleJsonArray();

            state._fsp--;

             current =iv_ruleJsonArray.getText(); 
            match(input,EOF,FOLLOW_2); 

            }

        }

            catch (RecognitionException re) {
                recover(input,re);
                appendSkippedTokens();
            }
        finally {
        }
        return current;
    }
    // $ANTLR end "entryRuleJsonArray"


    // $ANTLR start "ruleJsonArray"
    // InternalElkGraphJson.g:1857:1: ruleJsonArray returns [AntlrDatatypeRuleToken current=new AntlrDatatypeRuleToken()] : (kw= '[' (this_JsonValue_1= ruleJsonValue (kw= ',' this_JsonValue_3= ruleJsonValue )* )? (kw= ',' )? kw= ']' ) ;
    public final AntlrDatatypeRuleToken ruleJsonArray() throws RecognitionException {
        AntlrDatatypeRuleToken current = new AntlrDatatypeRuleToken();

        Token kw=null;
        AntlrDatatypeRuleToken this_JsonValue_1 = null;

        AntlrDatatypeRuleToken this_JsonValue_3 = null;



        	enterRule();

        try {
            // InternalElkGraphJson.g:1863:2: ( (kw= '[' (this_JsonValue_1= ruleJsonValue (kw= ',' this_JsonValue_3= ruleJsonValue )* )? (kw= ',' )? kw= ']' ) )
            // InternalElkGraphJson.g:1864:2: (kw= '[' (this_JsonValue_1= ruleJsonValue (kw= ',' this_JsonValue_3= ruleJsonValue )* )? (kw= ',' )? kw= ']' )
            {
            // InternalElkGraphJson.g:1864:2: (kw= '[' (this_JsonValue_1= ruleJsonValue (kw= ',' this_JsonValue_3= ruleJsonValue )* )? (kw= ',' )? kw= ']' )
            // InternalElkGraphJson.g:1865:3: kw= '[' (this_JsonValue_1= ruleJsonValue (kw= ',' this_JsonValue_3= ruleJsonValue )* )? (kw= ',' )? kw= ']'
            {
            kw=(Token)match(input,17,FOLLOW_21); 

            			current.merge(kw);
            			newLeafNode(kw, grammarAccess.getJsonArrayAccess().getLeftSquareBracketKeyword_0());
            		
            // InternalElkGraphJson.g:1870:3: (this_JsonValue_1= ruleJsonValue (kw= ',' this_JsonValue_3= ruleJsonValue )* )?
            int alt45=2;
            int LA45_0 = input.LA(1);

            if ( (LA45_0==RULE_STRING||(LA45_0>=RULE_SIGNED_INT && LA45_0<=RULE_FLOAT)||LA45_0==13||LA45_0==17||(LA45_0>=19 && LA45_0<=21)) ) {
                alt45=1;
            }
            switch (alt45) {
                case 1 :
                    // InternalElkGraphJson.g:1871:4: this_JsonValue_1= ruleJsonValue (kw= ',' this_JsonValue_3= ruleJsonValue )*
                    {

                    				newCompositeNode(grammarAccess.getJsonArrayAccess().getJsonValueParserRuleCall_1_0());
                    			
                    pushFollow(FOLLOW_14);
                    this_JsonValue_1=ruleJsonValue();

                    state._fsp--;


                    				current.merge(this_JsonValue_1);
                    			

                    				afterParserOrEnumRuleCall();
                    			
                    // InternalElkGraphJson.g:1881:4: (kw= ',' this_JsonValue_3= ruleJsonValue )*
                    loop44:
                    do {
                        int alt44=2;
                        int LA44_0 = input.LA(1);

                        if ( (LA44_0==14) ) {
                            int LA44_1 = input.LA(2);

                            if ( (LA44_1==RULE_STRING||(LA44_1>=RULE_SIGNED_INT && LA44_1<=RULE_FLOAT)||LA44_1==13||LA44_1==17||(LA44_1>=19 && LA44_1<=21)) ) {
                                alt44=1;
                            }


                        }


                        switch (alt44) {
                    	case 1 :
                    	    // InternalElkGraphJson.g:1882:5: kw= ',' this_JsonValue_3= ruleJsonValue
                    	    {
                    	    kw=(Token)match(input,14,FOLLOW_22); 

                    	    					current.merge(kw);
                    	    					newLeafNode(kw, grammarAccess.getJsonArrayAccess().getCommaKeyword_1_1_0());
                    	    				

                    	    					newCompositeNode(grammarAccess.getJsonArrayAccess().getJsonValueParserRuleCall_1_1_1());
                    	    				
                    	    pushFollow(FOLLOW_14);
                    	    this_JsonValue_3=ruleJsonValue();

                    	    state._fsp--;


                    	    					current.merge(this_JsonValue_3);
                    	    				

                    	    					afterParserOrEnumRuleCall();
                    	    				

                    	    }
                    	    break;

                    	default :
                    	    break loop44;
                        }
                    } while (true);


                    }
                    break;

            }

            // InternalElkGraphJson.g:1899:3: (kw= ',' )?
            int alt46=2;
            int LA46_0 = input.LA(1);

            if ( (LA46_0==14) ) {
                alt46=1;
            }
            switch (alt46) {
                case 1 :
                    // InternalElkGraphJson.g:1900:4: kw= ','
                    {
                    kw=(Token)match(input,14,FOLLOW_15); 

                    				current.merge(kw);
                    				newLeafNode(kw, grammarAccess.getJsonArrayAccess().getCommaKeyword_2());
                    			

                    }
                    break;

            }

            kw=(Token)match(input,18,FOLLOW_2); 

            			current.merge(kw);
            			newLeafNode(kw, grammarAccess.getJsonArrayAccess().getRightSquareBracketKeyword_3());
            		

            }


            }


            	leaveRule();

        }

            catch (RecognitionException re) {
                recover(input,re);
                appendSkippedTokens();
            }
        finally {
        }
        return current;
    }
    // $ANTLR end "ruleJsonArray"


    // $ANTLR start "entryRuleJsonMember"
    // InternalElkGraphJson.g:1915:1: entryRuleJsonMember returns [String current=null] : iv_ruleJsonMember= ruleJsonMember EOF ;
    public final String entryRuleJsonMember() throws RecognitionException {
        String current = null;

        AntlrDatatypeRuleToken iv_ruleJsonMember = null;


        try {
            // InternalElkGraphJson.g:1915:50: (iv_ruleJsonMember= ruleJsonMember EOF )
            // InternalElkGraphJson.g:1916:2: iv_ruleJsonMember= ruleJsonMember EOF
            {
             newCompositeNode(grammarAccess.getJsonMemberRule()); 
            pushFollow(FOLLOW_1);
            iv_ruleJsonMember=ruleJsonMember();

            state._fsp--;

             current =iv_ruleJsonMember.getText(); 
            match(input,EOF,FOLLOW_2); 

            }

        }

            catch (RecognitionException re) {
                recover(input,re);
                appendSkippedTokens();
            }
        finally {
        }
        return current;
    }
    // $ANTLR end "entryRuleJsonMember"


    // $ANTLR start "ruleJsonMember"
    // InternalElkGraphJson.g:1922:1: ruleJsonMember returns [AntlrDatatypeRuleToken current=new AntlrDatatypeRuleToken()] : ( (this_STRING_0= RULE_STRING | this_ID_1= RULE_ID ) kw= ':' this_JsonValue_3= ruleJsonValue ) ;
    public final AntlrDatatypeRuleToken ruleJsonMember() throws RecognitionException {
        AntlrDatatypeRuleToken current = new AntlrDatatypeRuleToken();

        Token this_STRING_0=null;
        Token this_ID_1=null;
        Token kw=null;
        AntlrDatatypeRuleToken this_JsonValue_3 = null;



        	enterRule();

        try {
            // InternalElkGraphJson.g:1928:2: ( ( (this_STRING_0= RULE_STRING | this_ID_1= RULE_ID ) kw= ':' this_JsonValue_3= ruleJsonValue ) )
            // InternalElkGraphJson.g:1929:2: ( (this_STRING_0= RULE_STRING | this_ID_1= RULE_ID ) kw= ':' this_JsonValue_3= ruleJsonValue )
            {
            // InternalElkGraphJson.g:1929:2: ( (this_STRING_0= RULE_STRING | this_ID_1= RULE_ID ) kw= ':' this_JsonValue_3= ruleJsonValue )
            // InternalElkGraphJson.g:1930:3: (this_STRING_0= RULE_STRING | this_ID_1= RULE_ID ) kw= ':' this_JsonValue_3= ruleJsonValue
            {
            // InternalElkGraphJson.g:1930:3: (this_STRING_0= RULE_STRING | this_ID_1= RULE_ID )
            int alt47=2;
            int LA47_0 = input.LA(1);

            if ( (LA47_0==RULE_STRING) ) {
                alt47=1;
            }
            else if ( (LA47_0==RULE_ID) ) {
                alt47=2;
            }
            else {
                NoViableAltException nvae =
                    new NoViableAltException("", 47, 0, input);

                throw nvae;
            }
            switch (alt47) {
                case 1 :
                    // InternalElkGraphJson.g:1931:4: this_STRING_0= RULE_STRING
                    {
                    this_STRING_0=(Token)match(input,RULE_STRING,FOLLOW_7); 

                    				current.merge(this_STRING_0);
                    			

                    				newLeafNode(this_STRING_0, grammarAccess.getJsonMemberAccess().getSTRINGTerminalRuleCall_0_0());
                    			

                    }
                    break;
                case 2 :
                    // InternalElkGraphJson.g:1939:4: this_ID_1= RULE_ID
                    {
                    this_ID_1=(Token)match(input,RULE_ID,FOLLOW_7); 

                    				current.merge(this_ID_1);
                    			

                    				newLeafNode(this_ID_1, grammarAccess.getJsonMemberAccess().getIDTerminalRuleCall_0_1());
                    			

                    }
                    break;

            }

            kw=(Token)match(input,16,FOLLOW_22); 

            			current.merge(kw);
            			newLeafNode(kw, grammarAccess.getJsonMemberAccess().getColonKeyword_1());
            		

            			newCompositeNode(grammarAccess.getJsonMemberAccess().getJsonValueParserRuleCall_2());
            		
            pushFollow(FOLLOW_2);
            this_JsonValue_3=ruleJsonValue();

            state._fsp--;


            			current.merge(this_JsonValue_3);
            		

            			afterParserOrEnumRuleCall();
            		

            }


            }


            	leaveRule();

        }

            catch (RecognitionException re) {
                recover(input,re);
                appendSkippedTokens();
            }
        finally {
        }
        return current;
    }
    // $ANTLR end "ruleJsonMember"


    // $ANTLR start "entryRuleJsonValue"
    // InternalElkGraphJson.g:1966:1: entryRuleJsonValue returns [String current=null] : iv_ruleJsonValue= ruleJsonValue EOF ;
    public final String entryRuleJsonValue() throws RecognitionException {
        String current = null;

        AntlrDatatypeRuleToken iv_ruleJsonValue = null;


        try {
            // InternalElkGraphJson.g:1966:49: (iv_ruleJsonValue= ruleJsonValue EOF )
            // InternalElkGraphJson.g:1967:2: iv_ruleJsonValue= ruleJsonValue EOF
            {
             newCompositeNode(grammarAccess.getJsonValueRule()); 
            pushFollow(FOLLOW_1);
            iv_ruleJsonValue=ruleJsonValue();

            state._fsp--;

             current =iv_ruleJsonValue.getText(); 
            match(input,EOF,FOLLOW_2); 

            }

        }

            catch (RecognitionException re) {
                recover(input,re);
                appendSkippedTokens();
            }
        finally {
        }
        return current;
    }
    // $ANTLR end "entryRuleJsonValue"


    // $ANTLR start "ruleJsonValue"
    // InternalElkGraphJson.g:1973:1: ruleJsonValue returns [AntlrDatatypeRuleToken current=new AntlrDatatypeRuleToken()] : (this_STRING_0= RULE_STRING | this_Number_1= ruleNumber | this_JsonObject_2= ruleJsonObject | this_JsonArray_3= ruleJsonArray | this_BooleanValue_4= ruleBooleanValue | kw= 'null' ) ;
    public final AntlrDatatypeRuleToken ruleJsonValue() throws RecognitionException {
        AntlrDatatypeRuleToken current = new AntlrDatatypeRuleToken();

        Token this_STRING_0=null;
        Token kw=null;
        AntlrDatatypeRuleToken this_Number_1 = null;

        AntlrDatatypeRuleToken this_JsonObject_2 = null;

        AntlrDatatypeRuleToken this_JsonArray_3 = null;

        AntlrDatatypeRuleToken this_BooleanValue_4 = null;



        	enterRule();

        try {
            // InternalElkGraphJson.g:1979:2: ( (this_STRING_0= RULE_STRING | this_Number_1= ruleNumber | this_JsonObject_2= ruleJsonObject | this_JsonArray_3= ruleJsonArray | this_BooleanValue_4= ruleBooleanValue | kw= 'null' ) )
            // InternalElkGraphJson.g:1980:2: (this_STRING_0= RULE_STRING | this_Number_1= ruleNumber | this_JsonObject_2= ruleJsonObject | this_JsonArray_3= ruleJsonArray | this_BooleanValue_4= ruleBooleanValue | kw= 'null' )
            {
            // InternalElkGraphJson.g:1980:2: (this_STRING_0= RULE_STRING | this_Number_1= ruleNumber | this_JsonObject_2= ruleJsonObject | this_JsonArray_3= ruleJsonArray | this_BooleanValue_4= ruleBooleanValue | kw= 'null' )
            int alt48=6;
            switch ( input.LA(1) ) {
            case RULE_STRING:
                {
                alt48=1;
                }
                break;
            case RULE_SIGNED_INT:
            case RULE_FLOAT:
                {
                alt48=2;
                }
                break;
            case 13:
                {
                alt48=3;
                }
                break;
            case 17:
                {
                alt48=4;
                }
                break;
            case 19:
            case 20:
                {
                alt48=5;
                }
                break;
            case 21:
                {
                alt48=6;
                }
                break;
            default:
                NoViableAltException nvae =
                    new NoViableAltException("", 48, 0, input);

                throw nvae;
            }

            switch (alt48) {
                case 1 :
                    // InternalElkGraphJson.g:1981:3: this_STRING_0= RULE_STRING
                    {
                    this_STRING_0=(Token)match(input,RULE_STRING,FOLLOW_2); 

                    			current.merge(this_STRING_0);
                    		

                    			newLeafNode(this_STRING_0, grammarAccess.getJsonValueAccess().getSTRINGTerminalRuleCall_0());
                    		

                    }
                    break;
                case 2 :
                    // InternalElkGraphJson.g:1989:3: this_Number_1= ruleNumber
                    {

                    			newCompositeNode(grammarAccess.getJsonValueAccess().getNumberParserRuleCall_1());
                    		
                    pushFollow(FOLLOW_2);
                    this_Number_1=ruleNumber();

                    state._fsp--;


                    			current.merge(this_Number_1);
                    		

                    			afterParserOrEnumRuleCall();
                    		

                    }
                    break;
                case 3 :
                    // InternalElkGraphJson.g:2000:3: this_JsonObject_2= ruleJsonObject
                    {

                    			newCompositeNode(grammarAccess.getJsonValueAccess().getJsonObjectParserRuleCall_2());
                    		
                    pushFollow(FOLLOW_2);
                    this_JsonObject_2=ruleJsonObject();

                    state._fsp--;


                    			current.merge(this_JsonObject_2);
                    		

                    			afterParserOrEnumRuleCall();
                    		

                    }
                    break;
                case 4 :
                    // InternalElkGraphJson.g:2011:3: this_JsonArray_3= ruleJsonArray
                    {

                    			newCompositeNode(grammarAccess.getJsonValueAccess().getJsonArrayParserRuleCall_3());
                    		
                    pushFollow(FOLLOW_2);
                    this_JsonArray_3=ruleJsonArray();

                    state._fsp--;


                    			current.merge(this_JsonArray_3);
                    		

                    			afterParserOrEnumRuleCall();
                    		

                    }
                    break;
                case 5 :
                    // InternalElkGraphJson.g:2022:3: this_BooleanValue_4= ruleBooleanValue
                    {

                    			newCompositeNode(grammarAccess.getJsonValueAccess().getBooleanValueParserRuleCall_4());
                    		
                    pushFollow(FOLLOW_2);
                    this_BooleanValue_4=ruleBooleanValue();

                    state._fsp--;


                    			current.merge(this_BooleanValue_4);
                    		

                    			afterParserOrEnumRuleCall();
                    		

                    }
                    break;
                case 6 :
                    // InternalElkGraphJson.g:2033:3: kw= 'null'
                    {
                    kw=(Token)match(input,21,FOLLOW_2); 

                    			current.merge(kw);
                    			newLeafNode(kw, grammarAccess.getJsonValueAccess().getNullKeyword_5());
                    		

                    }
                    break;

            }


            }


            	leaveRule();

        }

            catch (RecognitionException re) {
                recover(input,re);
                appendSkippedTokens();
            }
        finally {
        }
        return current;
    }
    // $ANTLR end "ruleJsonValue"


    // $ANTLR start "entryRuleKeyChildren"
    // InternalElkGraphJson.g:2042:1: entryRuleKeyChildren returns [String current=null] : iv_ruleKeyChildren= ruleKeyChildren EOF ;
    public final String entryRuleKeyChildren() throws RecognitionException {
        String current = null;

        AntlrDatatypeRuleToken iv_ruleKeyChildren = null;


        try {
            // InternalElkGraphJson.g:2042:51: (iv_ruleKeyChildren= ruleKeyChildren EOF )
            // InternalElkGraphJson.g:2043:2: iv_ruleKeyChildren= ruleKeyChildren EOF
            {
             newCompositeNode(grammarAccess.getKeyChildrenRule()); 
            pushFollow(FOLLOW_1);
            iv_ruleKeyChildren=ruleKeyChildren();

            state._fsp--;

             current =iv_ruleKeyChildren.getText(); 
            match(input,EOF,FOLLOW_2); 

            }

        }

            catch (RecognitionException re) {
                recover(input,re);
                appendSkippedTokens();
            }
        finally {
        }
        return current;
    }
    // $ANTLR end "entryRuleKeyChildren"


    // $ANTLR start "ruleKeyChildren"
    // InternalElkGraphJson.g:2049:1: ruleKeyChildren returns [AntlrDatatypeRuleToken current=new AntlrDatatypeRuleToken()] : (kw= '\"children\"' | kw= '\\'children\\'' | kw= 'children' ) ;
    public final AntlrDatatypeRuleToken ruleKeyChildren() throws RecognitionException {
        AntlrDatatypeRuleToken current = new AntlrDatatypeRuleToken();

        Token kw=null;


        	enterRule();

        try {
            // InternalElkGraphJson.g:2055:2: ( (kw= '\"children\"' | kw= '\\'children\\'' | kw= 'children' ) )
            // InternalElkGraphJson.g:2056:2: (kw= '\"children\"' | kw= '\\'children\\'' | kw= 'children' )
            {
            // InternalElkGraphJson.g:2056:2: (kw= '\"children\"' | kw= '\\'children\\'' | kw= 'children' )
            int alt49=3;
            switch ( input.LA(1) ) {
            case 22:
                {
                alt49=1;
                }
                break;
            case 23:
                {
                alt49=2;
                }
                break;
            case 24:
                {
                alt49=3;
                }
                break;
            default:
                NoViableAltException nvae =
                    new NoViableAltException("", 49, 0, input);

                throw nvae;
            }

            switch (alt49) {
                case 1 :
                    // InternalElkGraphJson.g:2057:3: kw= '\"children\"'
                    {
                    kw=(Token)match(input,22,FOLLOW_2); 

                    			current.merge(kw);
                    			newLeafNode(kw, grammarAccess.getKeyChildrenAccess().getChildrenKeyword_0());
                    		

                    }
                    break;
                case 2 :
                    // InternalElkGraphJson.g:2063:3: kw= '\\'children\\''
                    {
                    kw=(Token)match(input,23,FOLLOW_2); 

                    			current.merge(kw);
                    			newLeafNode(kw, grammarAccess.getKeyChildrenAccess().getChildrenKeyword_1());
                    		

                    }
                    break;
                case 3 :
                    // InternalElkGraphJson.g:2069:3: kw= 'children'
                    {
                    kw=(Token)match(input,24,FOLLOW_2); 

                    			current.merge(kw);
                    			newLeafNode(kw, grammarAccess.getKeyChildrenAccess().getChildrenKeyword_2());
                    		

                    }
                    break;

            }


            }


            	leaveRule();

        }

            catch (RecognitionException re) {
                recover(input,re);
                appendSkippedTokens();
            }
        finally {
        }
        return current;
    }
    // $ANTLR end "ruleKeyChildren"


    // $ANTLR start "entryRuleKeyPorts"
    // InternalElkGraphJson.g:2078:1: entryRuleKeyPorts returns [String current=null] : iv_ruleKeyPorts= ruleKeyPorts EOF ;
    public final String entryRuleKeyPorts() throws RecognitionException {
        String current = null;

        AntlrDatatypeRuleToken iv_ruleKeyPorts = null;


        try {
            // InternalElkGraphJson.g:2078:48: (iv_ruleKeyPorts= ruleKeyPorts EOF )
            // InternalElkGraphJson.g:2079:2: iv_ruleKeyPorts= ruleKeyPorts EOF
            {
             newCompositeNode(grammarAccess.getKeyPortsRule()); 
            pushFollow(FOLLOW_1);
            iv_ruleKeyPorts=ruleKeyPorts();

            state._fsp--;

             current =iv_ruleKeyPorts.getText(); 
            match(input,EOF,FOLLOW_2); 

            }

        }

            catch (RecognitionException re) {
                recover(input,re);
                appendSkippedTokens();
            }
        finally {
        }
        return current;
    }
    // $ANTLR end "entryRuleKeyPorts"


    // $ANTLR start "ruleKeyPorts"
    // InternalElkGraphJson.g:2085:1: ruleKeyPorts returns [AntlrDatatypeRuleToken current=new AntlrDatatypeRuleToken()] : (kw= '\"ports\"' | kw= '\\'ports\\'' | kw= 'ports' ) ;
    public final AntlrDatatypeRuleToken ruleKeyPorts() throws RecognitionException {
        AntlrDatatypeRuleToken current = new AntlrDatatypeRuleToken();

        Token kw=null;


        	enterRule();

        try {
            // InternalElkGraphJson.g:2091:2: ( (kw= '\"ports\"' | kw= '\\'ports\\'' | kw= 'ports' ) )
            // InternalElkGraphJson.g:2092:2: (kw= '\"ports\"' | kw= '\\'ports\\'' | kw= 'ports' )
            {
            // InternalElkGraphJson.g:2092:2: (kw= '\"ports\"' | kw= '\\'ports\\'' | kw= 'ports' )
            int alt50=3;
            switch ( input.LA(1) ) {
            case 25:
                {
                alt50=1;
                }
                break;
            case 26:
                {
                alt50=2;
                }
                break;
            case 27:
                {
                alt50=3;
                }
                break;
            default:
                NoViableAltException nvae =
                    new NoViableAltException("", 50, 0, input);

                throw nvae;
            }

            switch (alt50) {
                case 1 :
                    // InternalElkGraphJson.g:2093:3: kw= '\"ports\"'
                    {
                    kw=(Token)match(input,25,FOLLOW_2); 

                    			current.merge(kw);
                    			newLeafNode(kw, grammarAccess.getKeyPortsAccess().getPortsKeyword_0());
                    		

                    }
                    break;
                case 2 :
                    // InternalElkGraphJson.g:2099:3: kw= '\\'ports\\''
                    {
                    kw=(Token)match(input,26,FOLLOW_2); 

                    			current.merge(kw);
                    			newLeafNode(kw, grammarAccess.getKeyPortsAccess().getPortsKeyword_1());
                    		

                    }
                    break;
                case 3 :
                    // InternalElkGraphJson.g:2105:3: kw= 'ports'
                    {
                    kw=(Token)match(input,27,FOLLOW_2); 

                    			current.merge(kw);
                    			newLeafNode(kw, grammarAccess.getKeyPortsAccess().getPortsKeyword_2());
                    		

                    }
                    break;

            }


            }


            	leaveRule();

        }

            catch (RecognitionException re) {
                recover(input,re);
                appendSkippedTokens();
            }
        finally {
        }
        return current;
    }
    // $ANTLR end "ruleKeyPorts"


    // $ANTLR start "entryRuleKeyLabels"
    // InternalElkGraphJson.g:2114:1: entryRuleKeyLabels returns [String current=null] : iv_ruleKeyLabels= ruleKeyLabels EOF ;
    public final String entryRuleKeyLabels() throws RecognitionException {
        String current = null;

        AntlrDatatypeRuleToken iv_ruleKeyLabels = null;


        try {
            // InternalElkGraphJson.g:2114:49: (iv_ruleKeyLabels= ruleKeyLabels EOF )
            // InternalElkGraphJson.g:2115:2: iv_ruleKeyLabels= ruleKeyLabels EOF
            {
             newCompositeNode(grammarAccess.getKeyLabelsRule()); 
            pushFollow(FOLLOW_1);
            iv_ruleKeyLabels=ruleKeyLabels();

            state._fsp--;

             current =iv_ruleKeyLabels.getText(); 
            match(input,EOF,FOLLOW_2); 

            }

        }

            catch (RecognitionException re) {
                recover(input,re);
                appendSkippedTokens();
            }
        finally {
        }
        return current;
    }
    // $ANTLR end "entryRuleKeyLabels"


    // $ANTLR start "ruleKeyLabels"
    // InternalElkGraphJson.g:2121:1: ruleKeyLabels returns [AntlrDatatypeRuleToken current=new AntlrDatatypeRuleToken()] : (kw= '\"labels\"' | kw= '\\'labels\\'' | kw= 'labels' ) ;
    public final AntlrDatatypeRuleToken ruleKeyLabels() throws RecognitionException {
        AntlrDatatypeRuleToken current = new AntlrDatatypeRuleToken();

        Token kw=null;


        	enterRule();

        try {
            // InternalElkGraphJson.g:2127:2: ( (kw= '\"labels\"' | kw= '\\'labels\\'' | kw= 'labels' ) )
            // InternalElkGraphJson.g:2128:2: (kw= '\"labels\"' | kw= '\\'labels\\'' | kw= 'labels' )
            {
            // InternalElkGraphJson.g:2128:2: (kw= '\"labels\"' | kw= '\\'labels\\'' | kw= 'labels' )
            int alt51=3;
            switch ( input.LA(1) ) {
            case 28:
                {
                alt51=1;
                }
                break;
            case 29:
                {
                alt51=2;
                }
                break;
            case 30:
                {
                alt51=3;
                }
                break;
            default:
                NoViableAltException nvae =
                    new NoViableAltException("", 51, 0, input);

                throw nvae;
            }

            switch (alt51) {
                case 1 :
                    // InternalElkGraphJson.g:2129:3: kw= '\"labels\"'
                    {
                    kw=(Token)match(input,28,FOLLOW_2); 

                    			current.merge(kw);
                    			newLeafNode(kw, grammarAccess.getKeyLabelsAccess().getLabelsKeyword_0());
                    		

                    }
                    break;
                case 2 :
                    // InternalElkGraphJson.g:2135:3: kw= '\\'labels\\''
                    {
                    kw=(Token)match(input,29,FOLLOW_2); 

                    			current.merge(kw);
                    			newLeafNode(kw, grammarAccess.getKeyLabelsAccess().getLabelsKeyword_1());
                    		

                    }
                    break;
                case 3 :
                    // InternalElkGraphJson.g:2141:3: kw= 'labels'
                    {
                    kw=(Token)match(input,30,FOLLOW_2); 

                    			current.merge(kw);
                    			newLeafNode(kw, grammarAccess.getKeyLabelsAccess().getLabelsKeyword_2());
                    		

                    }
                    break;

            }


            }


            	leaveRule();

        }

            catch (RecognitionException re) {
                recover(input,re);
                appendSkippedTokens();
            }
        finally {
        }
        return current;
    }
    // $ANTLR end "ruleKeyLabels"


    // $ANTLR start "entryRuleKeyEdges"
    // InternalElkGraphJson.g:2150:1: entryRuleKeyEdges returns [String current=null] : iv_ruleKeyEdges= ruleKeyEdges EOF ;
    public final String entryRuleKeyEdges() throws RecognitionException {
        String current = null;

        AntlrDatatypeRuleToken iv_ruleKeyEdges = null;


        try {
            // InternalElkGraphJson.g:2150:48: (iv_ruleKeyEdges= ruleKeyEdges EOF )
            // InternalElkGraphJson.g:2151:2: iv_ruleKeyEdges= ruleKeyEdges EOF
            {
             newCompositeNode(grammarAccess.getKeyEdgesRule()); 
            pushFollow(FOLLOW_1);
            iv_ruleKeyEdges=ruleKeyEdges();

            state._fsp--;

             current =iv_ruleKeyEdges.getText(); 
            match(input,EOF,FOLLOW_2); 

            }

        }

            catch (RecognitionException re) {
                recover(input,re);
                appendSkippedTokens();
            }
        finally {
        }
        return current;
    }
    // $ANTLR end "entryRuleKeyEdges"


    // $ANTLR start "ruleKeyEdges"
    // InternalElkGraphJson.g:2157:1: ruleKeyEdges returns [AntlrDatatypeRuleToken current=new AntlrDatatypeRuleToken()] : (kw= '\"edges\"' | kw= '\\'edges\\'' | kw= 'edges' ) ;
    public final AntlrDatatypeRuleToken ruleKeyEdges() throws RecognitionException {
        AntlrDatatypeRuleToken current = new AntlrDatatypeRuleToken();

        Token kw=null;


        	enterRule();

        try {
            // InternalElkGraphJson.g:2163:2: ( (kw= '\"edges\"' | kw= '\\'edges\\'' | kw= 'edges' ) )
            // InternalElkGraphJson.g:2164:2: (kw= '\"edges\"' | kw= '\\'edges\\'' | kw= 'edges' )
            {
            // InternalElkGraphJson.g:2164:2: (kw= '\"edges\"' | kw= '\\'edges\\'' | kw= 'edges' )
            int alt52=3;
            switch ( input.LA(1) ) {
            case 31:
                {
                alt52=1;
                }
                break;
            case 32:
                {
                alt52=2;
                }
                break;
            case 33:
                {
                alt52=3;
                }
                break;
            default:
                NoViableAltException nvae =
                    new NoViableAltException("", 52, 0, input);

                throw nvae;
            }

            switch (alt52) {
                case 1 :
                    // InternalElkGraphJson.g:2165:3: kw= '\"edges\"'
                    {
                    kw=(Token)match(input,31,FOLLOW_2); 

                    			current.merge(kw);
                    			newLeafNode(kw, grammarAccess.getKeyEdgesAccess().getEdgesKeyword_0());
                    		

                    }
                    break;
                case 2 :
                    // InternalElkGraphJson.g:2171:3: kw= '\\'edges\\''
                    {
                    kw=(Token)match(input,32,FOLLOW_2); 

                    			current.merge(kw);
                    			newLeafNode(kw, grammarAccess.getKeyEdgesAccess().getEdgesKeyword_1());
                    		

                    }
                    break;
                case 3 :
                    // InternalElkGraphJson.g:2177:3: kw= 'edges'
                    {
                    kw=(Token)match(input,33,FOLLOW_2); 

                    			current.merge(kw);
                    			newLeafNode(kw, grammarAccess.getKeyEdgesAccess().getEdgesKeyword_2());
                    		

                    }
                    break;

            }


            }


            	leaveRule();

        }

            catch (RecognitionException re) {
                recover(input,re);
                appendSkippedTokens();
            }
        finally {
        }
        return current;
    }
    // $ANTLR end "ruleKeyEdges"


    // $ANTLR start "entryRuleKeyLayoutOptions"
    // InternalElkGraphJson.g:2186:1: entryRuleKeyLayoutOptions returns [String current=null] : iv_ruleKeyLayoutOptions= ruleKeyLayoutOptions EOF ;
    public final String entryRuleKeyLayoutOptions() throws RecognitionException {
        String current = null;

        AntlrDatatypeRuleToken iv_ruleKeyLayoutOptions = null;


        try {
            // InternalElkGraphJson.g:2186:56: (iv_ruleKeyLayoutOptions= ruleKeyLayoutOptions EOF )
            // InternalElkGraphJson.g:2187:2: iv_ruleKeyLayoutOptions= ruleKeyLayoutOptions EOF
            {
             newCompositeNode(grammarAccess.getKeyLayoutOptionsRule()); 
            pushFollow(FOLLOW_1);
            iv_ruleKeyLayoutOptions=ruleKeyLayoutOptions();

            state._fsp--;

             current =iv_ruleKeyLayoutOptions.getText(); 
            match(input,EOF,FOLLOW_2); 

            }

        }

            catch (RecognitionException re) {
                recover(input,re);
                appendSkippedTokens();
            }
        finally {
        }
        return current;
    }
    // $ANTLR end "entryRuleKeyLayoutOptions"


    // $ANTLR start "ruleKeyLayoutOptions"
    // InternalElkGraphJson.g:2193:1: ruleKeyLayoutOptions returns [AntlrDatatypeRuleToken current=new AntlrDatatypeRuleToken()] : (kw= '\"layoutOptions\"' | kw= '\\'layoutOptions\\'' | kw= 'layoutOptions' | kw= '\"properties\"' | kw= '\\'properties\\'' | kw= 'properties' ) ;
    public final AntlrDatatypeRuleToken ruleKeyLayoutOptions() throws RecognitionException {
        AntlrDatatypeRuleToken current = new AntlrDatatypeRuleToken();

        Token kw=null;


        	enterRule();

        try {
            // InternalElkGraphJson.g:2199:2: ( (kw= '\"layoutOptions\"' | kw= '\\'layoutOptions\\'' | kw= 'layoutOptions' | kw= '\"properties\"' | kw= '\\'properties\\'' | kw= 'properties' ) )
            // InternalElkGraphJson.g:2200:2: (kw= '\"layoutOptions\"' | kw= '\\'layoutOptions\\'' | kw= 'layoutOptions' | kw= '\"properties\"' | kw= '\\'properties\\'' | kw= 'properties' )
            {
            // InternalElkGraphJson.g:2200:2: (kw= '\"layoutOptions\"' | kw= '\\'layoutOptions\\'' | kw= 'layoutOptions' | kw= '\"properties\"' | kw= '\\'properties\\'' | kw= 'properties' )
            int alt53=6;
            switch ( input.LA(1) ) {
            case 34:
                {
                alt53=1;
                }
                break;
            case 35:
                {
                alt53=2;
                }
                break;
            case 36:
                {
                alt53=3;
                }
                break;
            case 37:
                {
                alt53=4;
                }
                break;
            case 38:
                {
                alt53=5;
                }
                break;
            case 39:
                {
                alt53=6;
                }
                break;
            default:
                NoViableAltException nvae =
                    new NoViableAltException("", 53, 0, input);

                throw nvae;
            }

            switch (alt53) {
                case 1 :
                    // InternalElkGraphJson.g:2201:3: kw= '\"layoutOptions\"'
                    {
                    kw=(Token)match(input,34,FOLLOW_2); 

                    			current.merge(kw);
                    			newLeafNode(kw, grammarAccess.getKeyLayoutOptionsAccess().getLayoutOptionsKeyword_0());
                    		

                    }
                    break;
                case 2 :
                    // InternalElkGraphJson.g:2207:3: kw= '\\'layoutOptions\\''
                    {
                    kw=(Token)match(input,35,FOLLOW_2); 

                    			current.merge(kw);
                    			newLeafNode(kw, grammarAccess.getKeyLayoutOptionsAccess().getLayoutOptionsKeyword_1());
                    		

                    }
                    break;
                case 3 :
                    // InternalElkGraphJson.g:2213:3: kw= 'layoutOptions'
                    {
                    kw=(Token)match(input,36,FOLLOW_2); 

                    			current.merge(kw);
                    			newLeafNode(kw, grammarAccess.getKeyLayoutOptionsAccess().getLayoutOptionsKeyword_2());
                    		

                    }
                    break;
                case 4 :
                    // InternalElkGraphJson.g:2219:3: kw= '\"properties\"'
                    {
                    kw=(Token)match(input,37,FOLLOW_2); 

                    			current.merge(kw);
                    			newLeafNode(kw, grammarAccess.getKeyLayoutOptionsAccess().getPropertiesKeyword_3());
                    		

                    }
                    break;
                case 5 :
                    // InternalElkGraphJson.g:2225:3: kw= '\\'properties\\''
                    {
                    kw=(Token)match(input,38,FOLLOW_2); 

                    			current.merge(kw);
                    			newLeafNode(kw, grammarAccess.getKeyLayoutOptionsAccess().getPropertiesKeyword_4());
                    		

                    }
                    break;
                case 6 :
                    // InternalElkGraphJson.g:2231:3: kw= 'properties'
                    {
                    kw=(Token)match(input,39,FOLLOW_2); 

                    			current.merge(kw);
                    			newLeafNode(kw, grammarAccess.getKeyLayoutOptionsAccess().getPropertiesKeyword_5());
                    		

                    }
                    break;

            }


            }


            	leaveRule();

        }

            catch (RecognitionException re) {
                recover(input,re);
                appendSkippedTokens();
            }
        finally {
        }
        return current;
    }
    // $ANTLR end "ruleKeyLayoutOptions"


    // $ANTLR start "entryRuleKeyId"
    // InternalElkGraphJson.g:2240:1: entryRuleKeyId returns [String current=null] : iv_ruleKeyId= ruleKeyId EOF ;
    public final String entryRuleKeyId() throws RecognitionException {
        String current = null;

        AntlrDatatypeRuleToken iv_ruleKeyId = null;


        try {
            // InternalElkGraphJson.g:2240:45: (iv_ruleKeyId= ruleKeyId EOF )
            // InternalElkGraphJson.g:2241:2: iv_ruleKeyId= ruleKeyId EOF
            {
             newCompositeNode(grammarAccess.getKeyIdRule()); 
            pushFollow(FOLLOW_1);
            iv_ruleKeyId=ruleKeyId();

            state._fsp--;

             current =iv_ruleKeyId.getText(); 
            match(input,EOF,FOLLOW_2); 

            }

        }

            catch (RecognitionException re) {
                recover(input,re);
                appendSkippedTokens();
            }
        finally {
        }
        return current;
    }
    // $ANTLR end "entryRuleKeyId"


    // $ANTLR start "ruleKeyId"
    // InternalElkGraphJson.g:2247:1: ruleKeyId returns [AntlrDatatypeRuleToken current=new AntlrDatatypeRuleToken()] : (kw= '\"id\"' | kw= '\\'id\\'' | kw= 'id' ) ;
    public final AntlrDatatypeRuleToken ruleKeyId() throws RecognitionException {
        AntlrDatatypeRuleToken current = new AntlrDatatypeRuleToken();

        Token kw=null;


        	enterRule();

        try {
            // InternalElkGraphJson.g:2253:2: ( (kw= '\"id\"' | kw= '\\'id\\'' | kw= 'id' ) )
            // InternalElkGraphJson.g:2254:2: (kw= '\"id\"' | kw= '\\'id\\'' | kw= 'id' )
            {
            // InternalElkGraphJson.g:2254:2: (kw= '\"id\"' | kw= '\\'id\\'' | kw= 'id' )
            int alt54=3;
            switch ( input.LA(1) ) {
            case 40:
                {
                alt54=1;
                }
                break;
            case 41:
                {
                alt54=2;
                }
                break;
            case 42:
                {
                alt54=3;
                }
                break;
            default:
                NoViableAltException nvae =
                    new NoViableAltException("", 54, 0, input);

                throw nvae;
            }

            switch (alt54) {
                case 1 :
                    // InternalElkGraphJson.g:2255:3: kw= '\"id\"'
                    {
                    kw=(Token)match(input,40,FOLLOW_2); 

                    			current.merge(kw);
                    			newLeafNode(kw, grammarAccess.getKeyIdAccess().getIdKeyword_0());
                    		

                    }
                    break;
                case 2 :
                    // InternalElkGraphJson.g:2261:3: kw= '\\'id\\''
                    {
                    kw=(Token)match(input,41,FOLLOW_2); 

                    			current.merge(kw);
                    			newLeafNode(kw, grammarAccess.getKeyIdAccess().getIdKeyword_1());
                    		

                    }
                    break;
                case 3 :
                    // InternalElkGraphJson.g:2267:3: kw= 'id'
                    {
                    kw=(Token)match(input,42,FOLLOW_2); 

                    			current.merge(kw);
                    			newLeafNode(kw, grammarAccess.getKeyIdAccess().getIdKeyword_2());
                    		

                    }
                    break;

            }


            }


            	leaveRule();

        }

            catch (RecognitionException re) {
                recover(input,re);
                appendSkippedTokens();
            }
        finally {
        }
        return current;
    }
    // $ANTLR end "ruleKeyId"


    // $ANTLR start "entryRuleKeyX"
    // InternalElkGraphJson.g:2276:1: entryRuleKeyX returns [String current=null] : iv_ruleKeyX= ruleKeyX EOF ;
    public final String entryRuleKeyX() throws RecognitionException {
        String current = null;

        AntlrDatatypeRuleToken iv_ruleKeyX = null;


        try {
            // InternalElkGraphJson.g:2276:44: (iv_ruleKeyX= ruleKeyX EOF )
            // InternalElkGraphJson.g:2277:2: iv_ruleKeyX= ruleKeyX EOF
            {
             newCompositeNode(grammarAccess.getKeyXRule()); 
            pushFollow(FOLLOW_1);
            iv_ruleKeyX=ruleKeyX();

            state._fsp--;

             current =iv_ruleKeyX.getText(); 
            match(input,EOF,FOLLOW_2); 

            }

        }

            catch (RecognitionException re) {
                recover(input,re);
                appendSkippedTokens();
            }
        finally {
        }
        return current;
    }
    // $ANTLR end "entryRuleKeyX"


    // $ANTLR start "ruleKeyX"
    // InternalElkGraphJson.g:2283:1: ruleKeyX returns [AntlrDatatypeRuleToken current=new AntlrDatatypeRuleToken()] : (kw= '\"x\"' | kw= '\\'x\\'' | kw= 'x' ) ;
    public final AntlrDatatypeRuleToken ruleKeyX() throws RecognitionException {
        AntlrDatatypeRuleToken current = new AntlrDatatypeRuleToken();

        Token kw=null;


        	enterRule();

        try {
            // InternalElkGraphJson.g:2289:2: ( (kw= '\"x\"' | kw= '\\'x\\'' | kw= 'x' ) )
            // InternalElkGraphJson.g:2290:2: (kw= '\"x\"' | kw= '\\'x\\'' | kw= 'x' )
            {
            // InternalElkGraphJson.g:2290:2: (kw= '\"x\"' | kw= '\\'x\\'' | kw= 'x' )
            int alt55=3;
            switch ( input.LA(1) ) {
            case 43:
                {
                alt55=1;
                }
                break;
            case 44:
                {
                alt55=2;
                }
                break;
            case 45:
                {
                alt55=3;
                }
                break;
            default:
                NoViableAltException nvae =
                    new NoViableAltException("", 55, 0, input);

                throw nvae;
            }

            switch (alt55) {
                case 1 :
                    // InternalElkGraphJson.g:2291:3: kw= '\"x\"'
                    {
                    kw=(Token)match(input,43,FOLLOW_2); 

                    			current.merge(kw);
                    			newLeafNode(kw, grammarAccess.getKeyXAccess().getXKeyword_0());
                    		

                    }
                    break;
                case 2 :
                    // InternalElkGraphJson.g:2297:3: kw= '\\'x\\''
                    {
                    kw=(Token)match(input,44,FOLLOW_2); 

                    			current.merge(kw);
                    			newLeafNode(kw, grammarAccess.getKeyXAccess().getXKeyword_1());
                    		

                    }
                    break;
                case 3 :
                    // InternalElkGraphJson.g:2303:3: kw= 'x'
                    {
                    kw=(Token)match(input,45,FOLLOW_2); 

                    			current.merge(kw);
                    			newLeafNode(kw, grammarAccess.getKeyXAccess().getXKeyword_2());
                    		

                    }
                    break;

            }


            }


            	leaveRule();

        }

            catch (RecognitionException re) {
                recover(input,re);
                appendSkippedTokens();
            }
        finally {
        }
        return current;
    }
    // $ANTLR end "ruleKeyX"


    // $ANTLR start "entryRuleKeyY"
    // InternalElkGraphJson.g:2312:1: entryRuleKeyY returns [String current=null] : iv_ruleKeyY= ruleKeyY EOF ;
    public final String entryRuleKeyY() throws RecognitionException {
        String current = null;

        AntlrDatatypeRuleToken iv_ruleKeyY = null;


        try {
            // InternalElkGraphJson.g:2312:44: (iv_ruleKeyY= ruleKeyY EOF )
            // InternalElkGraphJson.g:2313:2: iv_ruleKeyY= ruleKeyY EOF
            {
             newCompositeNode(grammarAccess.getKeyYRule()); 
            pushFollow(FOLLOW_1);
            iv_ruleKeyY=ruleKeyY();

            state._fsp--;

             current =iv_ruleKeyY.getText(); 
            match(input,EOF,FOLLOW_2); 

            }

        }

            catch (RecognitionException re) {
                recover(input,re);
                appendSkippedTokens();
            }
        finally {
        }
        return current;
    }
    // $ANTLR end "entryRuleKeyY"


    // $ANTLR start "ruleKeyY"
    // InternalElkGraphJson.g:2319:1: ruleKeyY returns [AntlrDatatypeRuleToken current=new AntlrDatatypeRuleToken()] : (kw= '\"y\"' | kw= '\\'y\\'' | kw= 'y' ) ;
    public final AntlrDatatypeRuleToken ruleKeyY() throws RecognitionException {
        AntlrDatatypeRuleToken current = new AntlrDatatypeRuleToken();

        Token kw=null;


        	enterRule();

        try {
            // InternalElkGraphJson.g:2325:2: ( (kw= '\"y\"' | kw= '\\'y\\'' | kw= 'y' ) )
            // InternalElkGraphJson.g:2326:2: (kw= '\"y\"' | kw= '\\'y\\'' | kw= 'y' )
            {
            // InternalElkGraphJson.g:2326:2: (kw= '\"y\"' | kw= '\\'y\\'' | kw= 'y' )
            int alt56=3;
            switch ( input.LA(1) ) {
            case 46:
                {
                alt56=1;
                }
                break;
            case 47:
                {
                alt56=2;
                }
                break;
            case 48:
                {
                alt56=3;
                }
                break;
            default:
                NoViableAltException nvae =
                    new NoViableAltException("", 56, 0, input);

                throw nvae;
            }

            switch (alt56) {
                case 1 :
                    // InternalElkGraphJson.g:2327:3: kw= '\"y\"'
                    {
                    kw=(Token)match(input,46,FOLLOW_2); 

                    			current.merge(kw);
                    			newLeafNode(kw, grammarAccess.getKeyYAccess().getYKeyword_0());
                    		

                    }
                    break;
                case 2 :
                    // InternalElkGraphJson.g:2333:3: kw= '\\'y\\''
                    {
                    kw=(Token)match(input,47,FOLLOW_2); 

                    			current.merge(kw);
                    			newLeafNode(kw, grammarAccess.getKeyYAccess().getYKeyword_1());
                    		

                    }
                    break;
                case 3 :
                    // InternalElkGraphJson.g:2339:3: kw= 'y'
                    {
                    kw=(Token)match(input,48,FOLLOW_2); 

                    			current.merge(kw);
                    			newLeafNode(kw, grammarAccess.getKeyYAccess().getYKeyword_2());
                    		

                    }
                    break;

            }


            }


            	leaveRule();

        }

            catch (RecognitionException re) {
                recover(input,re);
                appendSkippedTokens();
            }
        finally {
        }
        return current;
    }
    // $ANTLR end "ruleKeyY"


    // $ANTLR start "entryRuleKeyWidth"
    // InternalElkGraphJson.g:2348:1: entryRuleKeyWidth returns [String current=null] : iv_ruleKeyWidth= ruleKeyWidth EOF ;
    public final String entryRuleKeyWidth() throws RecognitionException {
        String current = null;

        AntlrDatatypeRuleToken iv_ruleKeyWidth = null;


        try {
            // InternalElkGraphJson.g:2348:48: (iv_ruleKeyWidth= ruleKeyWidth EOF )
            // InternalElkGraphJson.g:2349:2: iv_ruleKeyWidth= ruleKeyWidth EOF
            {
             newCompositeNode(grammarAccess.getKeyWidthRule()); 
            pushFollow(FOLLOW_1);
            iv_ruleKeyWidth=ruleKeyWidth();

            state._fsp--;

             current =iv_ruleKeyWidth.getText(); 
            match(input,EOF,FOLLOW_2); 

            }

        }

            catch (RecognitionException re) {
                recover(input,re);
                appendSkippedTokens();
            }
        finally {
        }
        return current;
    }
    // $ANTLR end "entryRuleKeyWidth"


    // $ANTLR start "ruleKeyWidth"
    // InternalElkGraphJson.g:2355:1: ruleKeyWidth returns [AntlrDatatypeRuleToken current=new AntlrDatatypeRuleToken()] : (kw= '\"width\"' | kw= '\\'width\\'' | kw= 'width' ) ;
    public final AntlrDatatypeRuleToken ruleKeyWidth() throws RecognitionException {
        AntlrDatatypeRuleToken current = new AntlrDatatypeRuleToken();

        Token kw=null;


        	enterRule();

        try {
            // InternalElkGraphJson.g:2361:2: ( (kw= '\"width\"' | kw= '\\'width\\'' | kw= 'width' ) )
            // InternalElkGraphJson.g:2362:2: (kw= '\"width\"' | kw= '\\'width\\'' | kw= 'width' )
            {
            // InternalElkGraphJson.g:2362:2: (kw= '\"width\"' | kw= '\\'width\\'' | kw= 'width' )
            int alt57=3;
            switch ( input.LA(1) ) {
            case 49:
                {
                alt57=1;
                }
                break;
            case 50:
                {
                alt57=2;
                }
                break;
            case 51:
                {
                alt57=3;
                }
                break;
            default:
                NoViableAltException nvae =
                    new NoViableAltException("", 57, 0, input);

                throw nvae;
            }

            switch (alt57) {
                case 1 :
                    // InternalElkGraphJson.g:2363:3: kw= '\"width\"'
                    {
                    kw=(Token)match(input,49,FOLLOW_2); 

                    			current.merge(kw);
                    			newLeafNode(kw, grammarAccess.getKeyWidthAccess().getWidthKeyword_0());
                    		

                    }
                    break;
                case 2 :
                    // InternalElkGraphJson.g:2369:3: kw= '\\'width\\''
                    {
                    kw=(Token)match(input,50,FOLLOW_2); 

                    			current.merge(kw);
                    			newLeafNode(kw, grammarAccess.getKeyWidthAccess().getWidthKeyword_1());
                    		

                    }
                    break;
                case 3 :
                    // InternalElkGraphJson.g:2375:3: kw= 'width'
                    {
                    kw=(Token)match(input,51,FOLLOW_2); 

                    			current.merge(kw);
                    			newLeafNode(kw, grammarAccess.getKeyWidthAccess().getWidthKeyword_2());
                    		

                    }
                    break;

            }


            }


            	leaveRule();

        }

            catch (RecognitionException re) {
                recover(input,re);
                appendSkippedTokens();
            }
        finally {
        }
        return current;
    }
    // $ANTLR end "ruleKeyWidth"


    // $ANTLR start "entryRuleKeyHeight"
    // InternalElkGraphJson.g:2384:1: entryRuleKeyHeight returns [String current=null] : iv_ruleKeyHeight= ruleKeyHeight EOF ;
    public final String entryRuleKeyHeight() throws RecognitionException {
        String current = null;

        AntlrDatatypeRuleToken iv_ruleKeyHeight = null;


        try {
            // InternalElkGraphJson.g:2384:49: (iv_ruleKeyHeight= ruleKeyHeight EOF )
            // InternalElkGraphJson.g:2385:2: iv_ruleKeyHeight= ruleKeyHeight EOF
            {
             newCompositeNode(grammarAccess.getKeyHeightRule()); 
            pushFollow(FOLLOW_1);
            iv_ruleKeyHeight=ruleKeyHeight();

            state._fsp--;

             current =iv_ruleKeyHeight.getText(); 
            match(input,EOF,FOLLOW_2); 

            }

        }

            catch (RecognitionException re) {
                recover(input,re);
                appendSkippedTokens();
            }
        finally {
        }
        return current;
    }
    // $ANTLR end "entryRuleKeyHeight"


    // $ANTLR start "ruleKeyHeight"
    // InternalElkGraphJson.g:2391:1: ruleKeyHeight returns [AntlrDatatypeRuleToken current=new AntlrDatatypeRuleToken()] : (kw= '\"height\"' | kw= '\\'height\\'' | kw= 'height' ) ;
    public final AntlrDatatypeRuleToken ruleKeyHeight() throws RecognitionException {
        AntlrDatatypeRuleToken current = new AntlrDatatypeRuleToken();

        Token kw=null;


        	enterRule();

        try {
            // InternalElkGraphJson.g:2397:2: ( (kw= '\"height\"' | kw= '\\'height\\'' | kw= 'height' ) )
            // InternalElkGraphJson.g:2398:2: (kw= '\"height\"' | kw= '\\'height\\'' | kw= 'height' )
            {
            // InternalElkGraphJson.g:2398:2: (kw= '\"height\"' | kw= '\\'height\\'' | kw= 'height' )
            int alt58=3;
            switch ( input.LA(1) ) {
            case 52:
                {
                alt58=1;
                }
                break;
            case 53:
                {
                alt58=2;
                }
                break;
            case 54:
                {
                alt58=3;
                }
                break;
            default:
                NoViableAltException nvae =
                    new NoViableAltException("", 58, 0, input);

                throw nvae;
            }

            switch (alt58) {
                case 1 :
                    // InternalElkGraphJson.g:2399:3: kw= '\"height\"'
                    {
                    kw=(Token)match(input,52,FOLLOW_2); 

                    			current.merge(kw);
                    			newLeafNode(kw, grammarAccess.getKeyHeightAccess().getHeightKeyword_0());
                    		

                    }
                    break;
                case 2 :
                    // InternalElkGraphJson.g:2405:3: kw= '\\'height\\''
                    {
                    kw=(Token)match(input,53,FOLLOW_2); 

                    			current.merge(kw);
                    			newLeafNode(kw, grammarAccess.getKeyHeightAccess().getHeightKeyword_1());
                    		

                    }
                    break;
                case 3 :
                    // InternalElkGraphJson.g:2411:3: kw= 'height'
                    {
                    kw=(Token)match(input,54,FOLLOW_2); 

                    			current.merge(kw);
                    			newLeafNode(kw, grammarAccess.getKeyHeightAccess().getHeightKeyword_2());
                    		

                    }
                    break;

            }


            }


            	leaveRule();

        }

            catch (RecognitionException re) {
                recover(input,re);
                appendSkippedTokens();
            }
        finally {
        }
        return current;
    }
    // $ANTLR end "ruleKeyHeight"


    // $ANTLR start "entryRuleKeySources"
    // InternalElkGraphJson.g:2420:1: entryRuleKeySources returns [String current=null] : iv_ruleKeySources= ruleKeySources EOF ;
    public final String entryRuleKeySources() throws RecognitionException {
        String current = null;

        AntlrDatatypeRuleToken iv_ruleKeySources = null;


        try {
            // InternalElkGraphJson.g:2420:50: (iv_ruleKeySources= ruleKeySources EOF )
            // InternalElkGraphJson.g:2421:2: iv_ruleKeySources= ruleKeySources EOF
            {
             newCompositeNode(grammarAccess.getKeySourcesRule()); 
            pushFollow(FOLLOW_1);
            iv_ruleKeySources=ruleKeySources();

            state._fsp--;

             current =iv_ruleKeySources.getText(); 
            match(input,EOF,FOLLOW_2); 

            }

        }

            catch (RecognitionException re) {
                recover(input,re);
                appendSkippedTokens();
            }
        finally {
        }
        return current;
    }
    // $ANTLR end "entryRuleKeySources"


    // $ANTLR start "ruleKeySources"
    // InternalElkGraphJson.g:2427:1: ruleKeySources returns [AntlrDatatypeRuleToken current=new AntlrDatatypeRuleToken()] : (kw= '\"sources\"' | kw= '\\'sources\\'' | kw= 'sources' ) ;
    public final AntlrDatatypeRuleToken ruleKeySources() throws RecognitionException {
        AntlrDatatypeRuleToken current = new AntlrDatatypeRuleToken();

        Token kw=null;


        	enterRule();

        try {
            // InternalElkGraphJson.g:2433:2: ( (kw= '\"sources\"' | kw= '\\'sources\\'' | kw= 'sources' ) )
            // InternalElkGraphJson.g:2434:2: (kw= '\"sources\"' | kw= '\\'sources\\'' | kw= 'sources' )
            {
            // InternalElkGraphJson.g:2434:2: (kw= '\"sources\"' | kw= '\\'sources\\'' | kw= 'sources' )
            int alt59=3;
            switch ( input.LA(1) ) {
            case 55:
                {
                alt59=1;
                }
                break;
            case 56:
                {
                alt59=2;
                }
                break;
            case 57:
                {
                alt59=3;
                }
                break;
            default:
                NoViableAltException nvae =
                    new NoViableAltException("", 59, 0, input);

                throw nvae;
            }

            switch (alt59) {
                case 1 :
                    // InternalElkGraphJson.g:2435:3: kw= '\"sources\"'
                    {
                    kw=(Token)match(input,55,FOLLOW_2); 

                    			current.merge(kw);
                    			newLeafNode(kw, grammarAccess.getKeySourcesAccess().getSourcesKeyword_0());
                    		

                    }
                    break;
                case 2 :
                    // InternalElkGraphJson.g:2441:3: kw= '\\'sources\\''
                    {
                    kw=(Token)match(input,56,FOLLOW_2); 

                    			current.merge(kw);
                    			newLeafNode(kw, grammarAccess.getKeySourcesAccess().getSourcesKeyword_1());
                    		

                    }
                    break;
                case 3 :
                    // InternalElkGraphJson.g:2447:3: kw= 'sources'
                    {
                    kw=(Token)match(input,57,FOLLOW_2); 

                    			current.merge(kw);
                    			newLeafNode(kw, grammarAccess.getKeySourcesAccess().getSourcesKeyword_2());
                    		

                    }
                    break;

            }


            }


            	leaveRule();

        }

            catch (RecognitionException re) {
                recover(input,re);
                appendSkippedTokens();
            }
        finally {
        }
        return current;
    }
    // $ANTLR end "ruleKeySources"


    // $ANTLR start "entryRuleKeyTargets"
    // InternalElkGraphJson.g:2456:1: entryRuleKeyTargets returns [String current=null] : iv_ruleKeyTargets= ruleKeyTargets EOF ;
    public final String entryRuleKeyTargets() throws RecognitionException {
        String current = null;

        AntlrDatatypeRuleToken iv_ruleKeyTargets = null;


        try {
            // InternalElkGraphJson.g:2456:50: (iv_ruleKeyTargets= ruleKeyTargets EOF )
            // InternalElkGraphJson.g:2457:2: iv_ruleKeyTargets= ruleKeyTargets EOF
            {
             newCompositeNode(grammarAccess.getKeyTargetsRule()); 
            pushFollow(FOLLOW_1);
            iv_ruleKeyTargets=ruleKeyTargets();

            state._fsp--;

             current =iv_ruleKeyTargets.getText(); 
            match(input,EOF,FOLLOW_2); 

            }

        }

            catch (RecognitionException re) {
                recover(input,re);
                appendSkippedTokens();
            }
        finally {
        }
        return current;
    }
    // $ANTLR end "entryRuleKeyTargets"


    // $ANTLR start "ruleKeyTargets"
    // InternalElkGraphJson.g:2463:1: ruleKeyTargets returns [AntlrDatatypeRuleToken current=new AntlrDatatypeRuleToken()] : (kw= '\"targets\"' | kw= '\\'targets\\'' | kw= 'targets' ) ;
    public final AntlrDatatypeRuleToken ruleKeyTargets() throws RecognitionException {
        AntlrDatatypeRuleToken current = new AntlrDatatypeRuleToken();

        Token kw=null;


        	enterRule();

        try {
            // InternalElkGraphJson.g:2469:2: ( (kw= '\"targets\"' | kw= '\\'targets\\'' | kw= 'targets' ) )
            // InternalElkGraphJson.g:2470:2: (kw= '\"targets\"' | kw= '\\'targets\\'' | kw= 'targets' )
            {
            // InternalElkGraphJson.g:2470:2: (kw= '\"targets\"' | kw= '\\'targets\\'' | kw= 'targets' )
            int alt60=3;
            switch ( input.LA(1) ) {
            case 58:
                {
                alt60=1;
                }
                break;
            case 59:
                {
                alt60=2;
                }
                break;
            case 60:
                {
                alt60=3;
                }
                break;
            default:
                NoViableAltException nvae =
                    new NoViableAltException("", 60, 0, input);

                throw nvae;
            }

            switch (alt60) {
                case 1 :
                    // InternalElkGraphJson.g:2471:3: kw= '\"targets\"'
                    {
                    kw=(Token)match(input,58,FOLLOW_2); 

                    			current.merge(kw);
                    			newLeafNode(kw, grammarAccess.getKeyTargetsAccess().getTargetsKeyword_0());
                    		

                    }
                    break;
                case 2 :
                    // InternalElkGraphJson.g:2477:3: kw= '\\'targets\\''
                    {
                    kw=(Token)match(input,59,FOLLOW_2); 

                    			current.merge(kw);
                    			newLeafNode(kw, grammarAccess.getKeyTargetsAccess().getTargetsKeyword_1());
                    		

                    }
                    break;
                case 3 :
                    // InternalElkGraphJson.g:2483:3: kw= 'targets'
                    {
                    kw=(Token)match(input,60,FOLLOW_2); 

                    			current.merge(kw);
                    			newLeafNode(kw, grammarAccess.getKeyTargetsAccess().getTargetsKeyword_2());
                    		

                    }
                    break;

            }


            }


            	leaveRule();

        }

            catch (RecognitionException re) {
                recover(input,re);
                appendSkippedTokens();
            }
        finally {
        }
        return current;
    }
    // $ANTLR end "ruleKeyTargets"


    // $ANTLR start "entryRuleKeyText"
    // InternalElkGraphJson.g:2492:1: entryRuleKeyText returns [String current=null] : iv_ruleKeyText= ruleKeyText EOF ;
    public final String entryRuleKeyText() throws RecognitionException {
        String current = null;

        AntlrDatatypeRuleToken iv_ruleKeyText = null;


        try {
            // InternalElkGraphJson.g:2492:47: (iv_ruleKeyText= ruleKeyText EOF )
            // InternalElkGraphJson.g:2493:2: iv_ruleKeyText= ruleKeyText EOF
            {
             newCompositeNode(grammarAccess.getKeyTextRule()); 
            pushFollow(FOLLOW_1);
            iv_ruleKeyText=ruleKeyText();

            state._fsp--;

             current =iv_ruleKeyText.getText(); 
            match(input,EOF,FOLLOW_2); 

            }

        }

            catch (RecognitionException re) {
                recover(input,re);
                appendSkippedTokens();
            }
        finally {
        }
        return current;
    }
    // $ANTLR end "entryRuleKeyText"


    // $ANTLR start "ruleKeyText"
    // InternalElkGraphJson.g:2499:1: ruleKeyText returns [AntlrDatatypeRuleToken current=new AntlrDatatypeRuleToken()] : (kw= '\"text\"' | kw= '\\'text\\'' | kw= 'text' ) ;
    public final AntlrDatatypeRuleToken ruleKeyText() throws RecognitionException {
        AntlrDatatypeRuleToken current = new AntlrDatatypeRuleToken();

        Token kw=null;


        	enterRule();

        try {
            // InternalElkGraphJson.g:2505:2: ( (kw= '\"text\"' | kw= '\\'text\\'' | kw= 'text' ) )
            // InternalElkGraphJson.g:2506:2: (kw= '\"text\"' | kw= '\\'text\\'' | kw= 'text' )
            {
            // InternalElkGraphJson.g:2506:2: (kw= '\"text\"' | kw= '\\'text\\'' | kw= 'text' )
            int alt61=3;
            switch ( input.LA(1) ) {
            case 61:
                {
                alt61=1;
                }
                break;
            case 62:
                {
                alt61=2;
                }
                break;
            case 63:
                {
                alt61=3;
                }
                break;
            default:
                NoViableAltException nvae =
                    new NoViableAltException("", 61, 0, input);

                throw nvae;
            }

            switch (alt61) {
                case 1 :
                    // InternalElkGraphJson.g:2507:3: kw= '\"text\"'
                    {
                    kw=(Token)match(input,61,FOLLOW_2); 

                    			current.merge(kw);
                    			newLeafNode(kw, grammarAccess.getKeyTextAccess().getTextKeyword_0());
                    		

                    }
                    break;
                case 2 :
                    // InternalElkGraphJson.g:2513:3: kw= '\\'text\\''
                    {
                    kw=(Token)match(input,62,FOLLOW_2); 

                    			current.merge(kw);
                    			newLeafNode(kw, grammarAccess.getKeyTextAccess().getTextKeyword_1());
                    		

                    }
                    break;
                case 3 :
                    // InternalElkGraphJson.g:2519:3: kw= 'text'
                    {
                    kw=(Token)match(input,63,FOLLOW_2); 

                    			current.merge(kw);
                    			newLeafNode(kw, grammarAccess.getKeyTextAccess().getTextKeyword_2());
                    		

                    }
                    break;

            }


            }


            	leaveRule();

        }

            catch (RecognitionException re) {
                recover(input,re);
                appendSkippedTokens();
            }
        finally {
        }
        return current;
    }
    // $ANTLR end "ruleKeyText"

    // Delegated rules


 

    public static final BitSet FOLLOW_1 = new BitSet(new long[]{0x0000000000000000L});
    public static final BitSet FOLLOW_2 = new BitSet(new long[]{0x0000000000000002L});
    public static final BitSet FOLLOW_3 = new BitSet(new long[]{0x007FFFFFFFC0C030L});
    public static final BitSet FOLLOW_4 = new BitSet(new long[]{0x000000000000C000L});
    public static final BitSet FOLLOW_5 = new BitSet(new long[]{0x007FFFFFFFC00030L});
    public static final BitSet FOLLOW_6 = new BitSet(new long[]{0x0000000000008000L});
    public static final BitSet FOLLOW_7 = new BitSet(new long[]{0x0000000000010000L});
    public static final BitSet FOLLOW_8 = new BitSet(new long[]{0x0000000000020000L});
    public static final BitSet FOLLOW_9 = new BitSet(new long[]{0x0000000000002000L});
    public static final BitSet FOLLOW_10 = new BitSet(new long[]{0xE07FFFFFFFC00030L});
    public static final BitSet FOLLOW_11 = new BitSet(new long[]{0x0000000000000010L});
    public static final BitSet FOLLOW_12 = new BitSet(new long[]{0x1FFFFFFFFFC00030L});
    public static final BitSet FOLLOW_13 = new BitSet(new long[]{0x0000000000044010L});
    public static final BitSet FOLLOW_14 = new BitSet(new long[]{0x0000000000044000L});
    public static final BitSet FOLLOW_15 = new BitSet(new long[]{0x0000000000040000L});
    public static final BitSet FOLLOW_16 = new BitSet(new long[]{0x0000000000046000L});
    public static final BitSet FOLLOW_17 = new BitSet(new long[]{0x000000000000C030L});
    public static final BitSet FOLLOW_18 = new BitSet(new long[]{0x0000000000000030L});
    public static final BitSet FOLLOW_19 = new BitSet(new long[]{0x00000000000000C0L});
    public static final BitSet FOLLOW_20 = new BitSet(new long[]{0x00000000001800D0L});
    public static final BitSet FOLLOW_21 = new BitSet(new long[]{0x00000000003E60D0L});
    public static final BitSet FOLLOW_22 = new BitSet(new long[]{0x00000000003A20D0L});

}