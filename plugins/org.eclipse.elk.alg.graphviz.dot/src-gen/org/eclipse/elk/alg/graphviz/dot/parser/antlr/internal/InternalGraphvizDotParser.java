/*******************************************************************************
 * Copyright (c) 2019 Kiel University and others.
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 * 
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.alg.graphviz.dot.parser.antlr.internal;

import org.eclipse.xtext.*;
import org.eclipse.xtext.parser.*;
import org.eclipse.xtext.parser.impl.*;
import org.eclipse.emf.ecore.util.EcoreUtil;
import org.eclipse.emf.ecore.EObject;
import org.eclipse.emf.common.util.Enumerator;
import org.eclipse.xtext.parser.antlr.AbstractInternalAntlrParser;
import org.eclipse.xtext.parser.antlr.XtextTokenStream;
import org.eclipse.xtext.parser.antlr.XtextTokenStream.HiddenTokens;
import org.eclipse.xtext.parser.antlr.AntlrDatatypeRuleToken;
import org.eclipse.elk.alg.graphviz.dot.services.GraphvizDotGrammarAccess;



import org.antlr.runtime.*;
import java.util.Stack;
import java.util.List;
import java.util.ArrayList;

/*******************************************************************************
 * Copyright (c) 2016 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
@SuppressWarnings("all")
public class InternalGraphvizDotParser extends AbstractInternalAntlrParser {
    public static final String[] tokenNames = new String[] {
        "<invalid>", "<EOR>", "<DOWN>", "<UP>", "RULE_ID", "RULE_INT", "RULE_FLOAT", "RULE_STRING", "RULE_PREC_LINE", "RULE_ML_COMMENT", "RULE_SL_COMMENT", "RULE_WS", "RULE_ANY_OTHER", "'strict'", "'{'", "'}'", "';'", "'='", "'['", "','", "']'", "':'", "'subgraph'", "'->'", "'--'", "'graph'", "'digraph'", "'node'", "'edge'"
    };
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


        public InternalGraphvizDotParser(TokenStream input) {
            this(input, new RecognizerSharedState());
        }
        public InternalGraphvizDotParser(TokenStream input, RecognizerSharedState state) {
            super(input, state);
             
        }
        

    public String[] getTokenNames() { return InternalGraphvizDotParser.tokenNames; }
    public String getGrammarFileName() { return "InternalGraphvizDot.g"; }



     	private GraphvizDotGrammarAccess grammarAccess;

        public InternalGraphvizDotParser(TokenStream input, GraphvizDotGrammarAccess grammarAccess) {
            this(input);
            this.grammarAccess = grammarAccess;
            registerRules(grammarAccess.getGrammar());
        }

        @Override
        protected String getFirstRuleName() {
        	return "GraphvizModel";
       	}

       	@Override
       	protected GraphvizDotGrammarAccess getGrammarAccess() {
       		return grammarAccess;
       	}




    // $ANTLR start "entryRuleGraphvizModel"
    // InternalGraphvizDot.g:72:1: entryRuleGraphvizModel returns [EObject current=null] : iv_ruleGraphvizModel= ruleGraphvizModel EOF ;
    public final EObject entryRuleGraphvizModel() throws RecognitionException {
        EObject current = null;

        EObject iv_ruleGraphvizModel = null;


        try {
            // InternalGraphvizDot.g:72:54: (iv_ruleGraphvizModel= ruleGraphvizModel EOF )
            // InternalGraphvizDot.g:73:2: iv_ruleGraphvizModel= ruleGraphvizModel EOF
            {
             newCompositeNode(grammarAccess.getGraphvizModelRule()); 
            pushFollow(FOLLOW_1);
            iv_ruleGraphvizModel=ruleGraphvizModel();

            state._fsp--;

             current =iv_ruleGraphvizModel; 
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
    // $ANTLR end "entryRuleGraphvizModel"


    // $ANTLR start "ruleGraphvizModel"
    // InternalGraphvizDot.g:79:1: ruleGraphvizModel returns [EObject current=null] : ( (lv_graphs_0_0= ruleGraph ) )* ;
    public final EObject ruleGraphvizModel() throws RecognitionException {
        EObject current = null;

        EObject lv_graphs_0_0 = null;



        	enterRule();

        try {
            // InternalGraphvizDot.g:85:2: ( ( (lv_graphs_0_0= ruleGraph ) )* )
            // InternalGraphvizDot.g:86:2: ( (lv_graphs_0_0= ruleGraph ) )*
            {
            // InternalGraphvizDot.g:86:2: ( (lv_graphs_0_0= ruleGraph ) )*
            loop1:
            do {
                int alt1=2;
                int LA1_0 = input.LA(1);

                if ( (LA1_0==13||(LA1_0>=25 && LA1_0<=26)) ) {
                    alt1=1;
                }


                switch (alt1) {
            	case 1 :
            	    // InternalGraphvizDot.g:87:3: (lv_graphs_0_0= ruleGraph )
            	    {
            	    // InternalGraphvizDot.g:87:3: (lv_graphs_0_0= ruleGraph )
            	    // InternalGraphvizDot.g:88:4: lv_graphs_0_0= ruleGraph
            	    {

            	    				newCompositeNode(grammarAccess.getGraphvizModelAccess().getGraphsGraphParserRuleCall_0());
            	    			
            	    pushFollow(FOLLOW_3);
            	    lv_graphs_0_0=ruleGraph();

            	    state._fsp--;


            	    				if (current==null) {
            	    					current = createModelElementForParent(grammarAccess.getGraphvizModelRule());
            	    				}
            	    				add(
            	    					current,
            	    					"graphs",
            	    					lv_graphs_0_0,
            	    					"org.eclipse.elk.alg.graphviz.dot.GraphvizDot.Graph");
            	    				afterParserOrEnumRuleCall();
            	    			

            	    }


            	    }
            	    break;

            	default :
            	    break loop1;
                }
            } while (true);


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
    // $ANTLR end "ruleGraphvizModel"


    // $ANTLR start "entryRuleGraph"
    // InternalGraphvizDot.g:108:1: entryRuleGraph returns [EObject current=null] : iv_ruleGraph= ruleGraph EOF ;
    public final EObject entryRuleGraph() throws RecognitionException {
        EObject current = null;

        EObject iv_ruleGraph = null;


        try {
            // InternalGraphvizDot.g:108:46: (iv_ruleGraph= ruleGraph EOF )
            // InternalGraphvizDot.g:109:2: iv_ruleGraph= ruleGraph EOF
            {
             newCompositeNode(grammarAccess.getGraphRule()); 
            pushFollow(FOLLOW_1);
            iv_ruleGraph=ruleGraph();

            state._fsp--;

             current =iv_ruleGraph; 
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
    // $ANTLR end "entryRuleGraph"


    // $ANTLR start "ruleGraph"
    // InternalGraphvizDot.g:115:1: ruleGraph returns [EObject current=null] : ( ( (lv_strict_0_0= 'strict' ) )? ( (lv_type_1_0= ruleGraphType ) ) ( (lv_name_2_0= ruleDotID ) )? otherlv_3= '{' ( (lv_statements_4_0= ruleStatement ) )* otherlv_5= '}' ) ;
    public final EObject ruleGraph() throws RecognitionException {
        EObject current = null;

        Token lv_strict_0_0=null;
        Token otherlv_3=null;
        Token otherlv_5=null;
        Enumerator lv_type_1_0 = null;

        AntlrDatatypeRuleToken lv_name_2_0 = null;

        EObject lv_statements_4_0 = null;



        	enterRule();

        try {
            // InternalGraphvizDot.g:121:2: ( ( ( (lv_strict_0_0= 'strict' ) )? ( (lv_type_1_0= ruleGraphType ) ) ( (lv_name_2_0= ruleDotID ) )? otherlv_3= '{' ( (lv_statements_4_0= ruleStatement ) )* otherlv_5= '}' ) )
            // InternalGraphvizDot.g:122:2: ( ( (lv_strict_0_0= 'strict' ) )? ( (lv_type_1_0= ruleGraphType ) ) ( (lv_name_2_0= ruleDotID ) )? otherlv_3= '{' ( (lv_statements_4_0= ruleStatement ) )* otherlv_5= '}' )
            {
            // InternalGraphvizDot.g:122:2: ( ( (lv_strict_0_0= 'strict' ) )? ( (lv_type_1_0= ruleGraphType ) ) ( (lv_name_2_0= ruleDotID ) )? otherlv_3= '{' ( (lv_statements_4_0= ruleStatement ) )* otherlv_5= '}' )
            // InternalGraphvizDot.g:123:3: ( (lv_strict_0_0= 'strict' ) )? ( (lv_type_1_0= ruleGraphType ) ) ( (lv_name_2_0= ruleDotID ) )? otherlv_3= '{' ( (lv_statements_4_0= ruleStatement ) )* otherlv_5= '}'
            {
            // InternalGraphvizDot.g:123:3: ( (lv_strict_0_0= 'strict' ) )?
            int alt2=2;
            int LA2_0 = input.LA(1);

            if ( (LA2_0==13) ) {
                alt2=1;
            }
            switch (alt2) {
                case 1 :
                    // InternalGraphvizDot.g:124:4: (lv_strict_0_0= 'strict' )
                    {
                    // InternalGraphvizDot.g:124:4: (lv_strict_0_0= 'strict' )
                    // InternalGraphvizDot.g:125:5: lv_strict_0_0= 'strict'
                    {
                    lv_strict_0_0=(Token)match(input,13,FOLLOW_4); 

                    					newLeafNode(lv_strict_0_0, grammarAccess.getGraphAccess().getStrictStrictKeyword_0_0());
                    				

                    					if (current==null) {
                    						current = createModelElement(grammarAccess.getGraphRule());
                    					}
                    					setWithLastConsumed(current, "strict", true, "strict");
                    				

                    }


                    }
                    break;

            }

            // InternalGraphvizDot.g:137:3: ( (lv_type_1_0= ruleGraphType ) )
            // InternalGraphvizDot.g:138:4: (lv_type_1_0= ruleGraphType )
            {
            // InternalGraphvizDot.g:138:4: (lv_type_1_0= ruleGraphType )
            // InternalGraphvizDot.g:139:5: lv_type_1_0= ruleGraphType
            {

            					newCompositeNode(grammarAccess.getGraphAccess().getTypeGraphTypeEnumRuleCall_1_0());
            				
            pushFollow(FOLLOW_5);
            lv_type_1_0=ruleGraphType();

            state._fsp--;


            					if (current==null) {
            						current = createModelElementForParent(grammarAccess.getGraphRule());
            					}
            					set(
            						current,
            						"type",
            						lv_type_1_0,
            						"org.eclipse.elk.alg.graphviz.dot.GraphvizDot.GraphType");
            					afterParserOrEnumRuleCall();
            				

            }


            }

            // InternalGraphvizDot.g:156:3: ( (lv_name_2_0= ruleDotID ) )?
            int alt3=2;
            int LA3_0 = input.LA(1);

            if ( ((LA3_0>=RULE_ID && LA3_0<=RULE_STRING)) ) {
                alt3=1;
            }
            switch (alt3) {
                case 1 :
                    // InternalGraphvizDot.g:157:4: (lv_name_2_0= ruleDotID )
                    {
                    // InternalGraphvizDot.g:157:4: (lv_name_2_0= ruleDotID )
                    // InternalGraphvizDot.g:158:5: lv_name_2_0= ruleDotID
                    {

                    					newCompositeNode(grammarAccess.getGraphAccess().getNameDotIDParserRuleCall_2_0());
                    				
                    pushFollow(FOLLOW_6);
                    lv_name_2_0=ruleDotID();

                    state._fsp--;


                    					if (current==null) {
                    						current = createModelElementForParent(grammarAccess.getGraphRule());
                    					}
                    					set(
                    						current,
                    						"name",
                    						lv_name_2_0,
                    						"org.eclipse.elk.alg.graphviz.dot.GraphvizDot.DotID");
                    					afterParserOrEnumRuleCall();
                    				

                    }


                    }
                    break;

            }

            otherlv_3=(Token)match(input,14,FOLLOW_7); 

            			newLeafNode(otherlv_3, grammarAccess.getGraphAccess().getLeftCurlyBracketKeyword_3());
            		
            // InternalGraphvizDot.g:179:3: ( (lv_statements_4_0= ruleStatement ) )*
            loop4:
            do {
                int alt4=2;
                int LA4_0 = input.LA(1);

                if ( ((LA4_0>=RULE_ID && LA4_0<=RULE_STRING)||LA4_0==14||LA4_0==22||LA4_0==25||(LA4_0>=27 && LA4_0<=28)) ) {
                    alt4=1;
                }


                switch (alt4) {
            	case 1 :
            	    // InternalGraphvizDot.g:180:4: (lv_statements_4_0= ruleStatement )
            	    {
            	    // InternalGraphvizDot.g:180:4: (lv_statements_4_0= ruleStatement )
            	    // InternalGraphvizDot.g:181:5: lv_statements_4_0= ruleStatement
            	    {

            	    					newCompositeNode(grammarAccess.getGraphAccess().getStatementsStatementParserRuleCall_4_0());
            	    				
            	    pushFollow(FOLLOW_7);
            	    lv_statements_4_0=ruleStatement();

            	    state._fsp--;


            	    					if (current==null) {
            	    						current = createModelElementForParent(grammarAccess.getGraphRule());
            	    					}
            	    					add(
            	    						current,
            	    						"statements",
            	    						lv_statements_4_0,
            	    						"org.eclipse.elk.alg.graphviz.dot.GraphvizDot.Statement");
            	    					afterParserOrEnumRuleCall();
            	    				

            	    }


            	    }
            	    break;

            	default :
            	    break loop4;
                }
            } while (true);

            otherlv_5=(Token)match(input,15,FOLLOW_2); 

            			newLeafNode(otherlv_5, grammarAccess.getGraphAccess().getRightCurlyBracketKeyword_5());
            		

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
    // $ANTLR end "ruleGraph"


    // $ANTLR start "entryRuleStatement"
    // InternalGraphvizDot.g:206:1: entryRuleStatement returns [EObject current=null] : iv_ruleStatement= ruleStatement EOF ;
    public final EObject entryRuleStatement() throws RecognitionException {
        EObject current = null;

        EObject iv_ruleStatement = null;


        try {
            // InternalGraphvizDot.g:206:50: (iv_ruleStatement= ruleStatement EOF )
            // InternalGraphvizDot.g:207:2: iv_ruleStatement= ruleStatement EOF
            {
             newCompositeNode(grammarAccess.getStatementRule()); 
            pushFollow(FOLLOW_1);
            iv_ruleStatement=ruleStatement();

            state._fsp--;

             current =iv_ruleStatement; 
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
    // $ANTLR end "entryRuleStatement"


    // $ANTLR start "ruleStatement"
    // InternalGraphvizDot.g:213:1: ruleStatement returns [EObject current=null] : ( (this_NodeStatement_0= ruleNodeStatement | this_EdgeStatement_1= ruleEdgeStatement | this_AttributeStatement_2= ruleAttributeStatement | this_Attribute_3= ruleAttribute | this_Subgraph_4= ruleSubgraph ) (otherlv_5= ';' )? ) ;
    public final EObject ruleStatement() throws RecognitionException {
        EObject current = null;

        Token otherlv_5=null;
        EObject this_NodeStatement_0 = null;

        EObject this_EdgeStatement_1 = null;

        EObject this_AttributeStatement_2 = null;

        EObject this_Attribute_3 = null;

        EObject this_Subgraph_4 = null;



        	enterRule();

        try {
            // InternalGraphvizDot.g:219:2: ( ( (this_NodeStatement_0= ruleNodeStatement | this_EdgeStatement_1= ruleEdgeStatement | this_AttributeStatement_2= ruleAttributeStatement | this_Attribute_3= ruleAttribute | this_Subgraph_4= ruleSubgraph ) (otherlv_5= ';' )? ) )
            // InternalGraphvizDot.g:220:2: ( (this_NodeStatement_0= ruleNodeStatement | this_EdgeStatement_1= ruleEdgeStatement | this_AttributeStatement_2= ruleAttributeStatement | this_Attribute_3= ruleAttribute | this_Subgraph_4= ruleSubgraph ) (otherlv_5= ';' )? )
            {
            // InternalGraphvizDot.g:220:2: ( (this_NodeStatement_0= ruleNodeStatement | this_EdgeStatement_1= ruleEdgeStatement | this_AttributeStatement_2= ruleAttributeStatement | this_Attribute_3= ruleAttribute | this_Subgraph_4= ruleSubgraph ) (otherlv_5= ';' )? )
            // InternalGraphvizDot.g:221:3: (this_NodeStatement_0= ruleNodeStatement | this_EdgeStatement_1= ruleEdgeStatement | this_AttributeStatement_2= ruleAttributeStatement | this_Attribute_3= ruleAttribute | this_Subgraph_4= ruleSubgraph ) (otherlv_5= ';' )?
            {
            // InternalGraphvizDot.g:221:3: (this_NodeStatement_0= ruleNodeStatement | this_EdgeStatement_1= ruleEdgeStatement | this_AttributeStatement_2= ruleAttributeStatement | this_Attribute_3= ruleAttribute | this_Subgraph_4= ruleSubgraph )
            int alt5=5;
            alt5 = dfa5.predict(input);
            switch (alt5) {
                case 1 :
                    // InternalGraphvizDot.g:222:4: this_NodeStatement_0= ruleNodeStatement
                    {

                    				newCompositeNode(grammarAccess.getStatementAccess().getNodeStatementParserRuleCall_0_0());
                    			
                    pushFollow(FOLLOW_8);
                    this_NodeStatement_0=ruleNodeStatement();

                    state._fsp--;


                    				current = this_NodeStatement_0;
                    				afterParserOrEnumRuleCall();
                    			

                    }
                    break;
                case 2 :
                    // InternalGraphvizDot.g:231:4: this_EdgeStatement_1= ruleEdgeStatement
                    {

                    				newCompositeNode(grammarAccess.getStatementAccess().getEdgeStatementParserRuleCall_0_1());
                    			
                    pushFollow(FOLLOW_8);
                    this_EdgeStatement_1=ruleEdgeStatement();

                    state._fsp--;


                    				current = this_EdgeStatement_1;
                    				afterParserOrEnumRuleCall();
                    			

                    }
                    break;
                case 3 :
                    // InternalGraphvizDot.g:240:4: this_AttributeStatement_2= ruleAttributeStatement
                    {

                    				newCompositeNode(grammarAccess.getStatementAccess().getAttributeStatementParserRuleCall_0_2());
                    			
                    pushFollow(FOLLOW_8);
                    this_AttributeStatement_2=ruleAttributeStatement();

                    state._fsp--;


                    				current = this_AttributeStatement_2;
                    				afterParserOrEnumRuleCall();
                    			

                    }
                    break;
                case 4 :
                    // InternalGraphvizDot.g:249:4: this_Attribute_3= ruleAttribute
                    {

                    				newCompositeNode(grammarAccess.getStatementAccess().getAttributeParserRuleCall_0_3());
                    			
                    pushFollow(FOLLOW_8);
                    this_Attribute_3=ruleAttribute();

                    state._fsp--;


                    				current = this_Attribute_3;
                    				afterParserOrEnumRuleCall();
                    			

                    }
                    break;
                case 5 :
                    // InternalGraphvizDot.g:258:4: this_Subgraph_4= ruleSubgraph
                    {

                    				newCompositeNode(grammarAccess.getStatementAccess().getSubgraphParserRuleCall_0_4());
                    			
                    pushFollow(FOLLOW_8);
                    this_Subgraph_4=ruleSubgraph();

                    state._fsp--;


                    				current = this_Subgraph_4;
                    				afterParserOrEnumRuleCall();
                    			

                    }
                    break;

            }

            // InternalGraphvizDot.g:267:3: (otherlv_5= ';' )?
            int alt6=2;
            int LA6_0 = input.LA(1);

            if ( (LA6_0==16) ) {
                alt6=1;
            }
            switch (alt6) {
                case 1 :
                    // InternalGraphvizDot.g:268:4: otherlv_5= ';'
                    {
                    otherlv_5=(Token)match(input,16,FOLLOW_2); 

                    				newLeafNode(otherlv_5, grammarAccess.getStatementAccess().getSemicolonKeyword_1());
                    			

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
    // $ANTLR end "ruleStatement"


    // $ANTLR start "entryRuleAttribute"
    // InternalGraphvizDot.g:277:1: entryRuleAttribute returns [EObject current=null] : iv_ruleAttribute= ruleAttribute EOF ;
    public final EObject entryRuleAttribute() throws RecognitionException {
        EObject current = null;

        EObject iv_ruleAttribute = null;


        try {
            // InternalGraphvizDot.g:277:50: (iv_ruleAttribute= ruleAttribute EOF )
            // InternalGraphvizDot.g:278:2: iv_ruleAttribute= ruleAttribute EOF
            {
             newCompositeNode(grammarAccess.getAttributeRule()); 
            pushFollow(FOLLOW_1);
            iv_ruleAttribute=ruleAttribute();

            state._fsp--;

             current =iv_ruleAttribute; 
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
    // $ANTLR end "entryRuleAttribute"


    // $ANTLR start "ruleAttribute"
    // InternalGraphvizDot.g:284:1: ruleAttribute returns [EObject current=null] : ( ( (lv_name_0_0= ruleDotID ) ) otherlv_1= '=' ( (lv_value_2_0= ruleDotID ) ) ) ;
    public final EObject ruleAttribute() throws RecognitionException {
        EObject current = null;

        Token otherlv_1=null;
        AntlrDatatypeRuleToken lv_name_0_0 = null;

        AntlrDatatypeRuleToken lv_value_2_0 = null;



        	enterRule();

        try {
            // InternalGraphvizDot.g:290:2: ( ( ( (lv_name_0_0= ruleDotID ) ) otherlv_1= '=' ( (lv_value_2_0= ruleDotID ) ) ) )
            // InternalGraphvizDot.g:291:2: ( ( (lv_name_0_0= ruleDotID ) ) otherlv_1= '=' ( (lv_value_2_0= ruleDotID ) ) )
            {
            // InternalGraphvizDot.g:291:2: ( ( (lv_name_0_0= ruleDotID ) ) otherlv_1= '=' ( (lv_value_2_0= ruleDotID ) ) )
            // InternalGraphvizDot.g:292:3: ( (lv_name_0_0= ruleDotID ) ) otherlv_1= '=' ( (lv_value_2_0= ruleDotID ) )
            {
            // InternalGraphvizDot.g:292:3: ( (lv_name_0_0= ruleDotID ) )
            // InternalGraphvizDot.g:293:4: (lv_name_0_0= ruleDotID )
            {
            // InternalGraphvizDot.g:293:4: (lv_name_0_0= ruleDotID )
            // InternalGraphvizDot.g:294:5: lv_name_0_0= ruleDotID
            {

            					newCompositeNode(grammarAccess.getAttributeAccess().getNameDotIDParserRuleCall_0_0());
            				
            pushFollow(FOLLOW_9);
            lv_name_0_0=ruleDotID();

            state._fsp--;


            					if (current==null) {
            						current = createModelElementForParent(grammarAccess.getAttributeRule());
            					}
            					set(
            						current,
            						"name",
            						lv_name_0_0,
            						"org.eclipse.elk.alg.graphviz.dot.GraphvizDot.DotID");
            					afterParserOrEnumRuleCall();
            				

            }


            }

            otherlv_1=(Token)match(input,17,FOLLOW_10); 

            			newLeafNode(otherlv_1, grammarAccess.getAttributeAccess().getEqualsSignKeyword_1());
            		
            // InternalGraphvizDot.g:315:3: ( (lv_value_2_0= ruleDotID ) )
            // InternalGraphvizDot.g:316:4: (lv_value_2_0= ruleDotID )
            {
            // InternalGraphvizDot.g:316:4: (lv_value_2_0= ruleDotID )
            // InternalGraphvizDot.g:317:5: lv_value_2_0= ruleDotID
            {

            					newCompositeNode(grammarAccess.getAttributeAccess().getValueDotIDParserRuleCall_2_0());
            				
            pushFollow(FOLLOW_2);
            lv_value_2_0=ruleDotID();

            state._fsp--;


            					if (current==null) {
            						current = createModelElementForParent(grammarAccess.getAttributeRule());
            					}
            					set(
            						current,
            						"value",
            						lv_value_2_0,
            						"org.eclipse.elk.alg.graphviz.dot.GraphvizDot.DotID");
            					afterParserOrEnumRuleCall();
            				

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
    // $ANTLR end "ruleAttribute"


    // $ANTLR start "entryRuleNodeStatement"
    // InternalGraphvizDot.g:338:1: entryRuleNodeStatement returns [EObject current=null] : iv_ruleNodeStatement= ruleNodeStatement EOF ;
    public final EObject entryRuleNodeStatement() throws RecognitionException {
        EObject current = null;

        EObject iv_ruleNodeStatement = null;


        try {
            // InternalGraphvizDot.g:338:54: (iv_ruleNodeStatement= ruleNodeStatement EOF )
            // InternalGraphvizDot.g:339:2: iv_ruleNodeStatement= ruleNodeStatement EOF
            {
             newCompositeNode(grammarAccess.getNodeStatementRule()); 
            pushFollow(FOLLOW_1);
            iv_ruleNodeStatement=ruleNodeStatement();

            state._fsp--;

             current =iv_ruleNodeStatement; 
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
    // $ANTLR end "entryRuleNodeStatement"


    // $ANTLR start "ruleNodeStatement"
    // InternalGraphvizDot.g:345:1: ruleNodeStatement returns [EObject current=null] : ( ( (lv_node_0_0= ruleNode ) ) (otherlv_1= '[' ( ( (lv_attributes_2_0= ruleListAttribute ) ) ( (otherlv_3= ',' )? ( (lv_attributes_4_0= ruleListAttribute ) ) )* )? otherlv_5= ']' )? ) ;
    public final EObject ruleNodeStatement() throws RecognitionException {
        EObject current = null;

        Token otherlv_1=null;
        Token otherlv_3=null;
        Token otherlv_5=null;
        EObject lv_node_0_0 = null;

        EObject lv_attributes_2_0 = null;

        EObject lv_attributes_4_0 = null;



        	enterRule();

        try {
            // InternalGraphvizDot.g:351:2: ( ( ( (lv_node_0_0= ruleNode ) ) (otherlv_1= '[' ( ( (lv_attributes_2_0= ruleListAttribute ) ) ( (otherlv_3= ',' )? ( (lv_attributes_4_0= ruleListAttribute ) ) )* )? otherlv_5= ']' )? ) )
            // InternalGraphvizDot.g:352:2: ( ( (lv_node_0_0= ruleNode ) ) (otherlv_1= '[' ( ( (lv_attributes_2_0= ruleListAttribute ) ) ( (otherlv_3= ',' )? ( (lv_attributes_4_0= ruleListAttribute ) ) )* )? otherlv_5= ']' )? )
            {
            // InternalGraphvizDot.g:352:2: ( ( (lv_node_0_0= ruleNode ) ) (otherlv_1= '[' ( ( (lv_attributes_2_0= ruleListAttribute ) ) ( (otherlv_3= ',' )? ( (lv_attributes_4_0= ruleListAttribute ) ) )* )? otherlv_5= ']' )? )
            // InternalGraphvizDot.g:353:3: ( (lv_node_0_0= ruleNode ) ) (otherlv_1= '[' ( ( (lv_attributes_2_0= ruleListAttribute ) ) ( (otherlv_3= ',' )? ( (lv_attributes_4_0= ruleListAttribute ) ) )* )? otherlv_5= ']' )?
            {
            // InternalGraphvizDot.g:353:3: ( (lv_node_0_0= ruleNode ) )
            // InternalGraphvizDot.g:354:4: (lv_node_0_0= ruleNode )
            {
            // InternalGraphvizDot.g:354:4: (lv_node_0_0= ruleNode )
            // InternalGraphvizDot.g:355:5: lv_node_0_0= ruleNode
            {

            					newCompositeNode(grammarAccess.getNodeStatementAccess().getNodeNodeParserRuleCall_0_0());
            				
            pushFollow(FOLLOW_11);
            lv_node_0_0=ruleNode();

            state._fsp--;


            					if (current==null) {
            						current = createModelElementForParent(grammarAccess.getNodeStatementRule());
            					}
            					set(
            						current,
            						"node",
            						lv_node_0_0,
            						"org.eclipse.elk.alg.graphviz.dot.GraphvizDot.Node");
            					afterParserOrEnumRuleCall();
            				

            }


            }

            // InternalGraphvizDot.g:372:3: (otherlv_1= '[' ( ( (lv_attributes_2_0= ruleListAttribute ) ) ( (otherlv_3= ',' )? ( (lv_attributes_4_0= ruleListAttribute ) ) )* )? otherlv_5= ']' )?
            int alt10=2;
            int LA10_0 = input.LA(1);

            if ( (LA10_0==18) ) {
                alt10=1;
            }
            switch (alt10) {
                case 1 :
                    // InternalGraphvizDot.g:373:4: otherlv_1= '[' ( ( (lv_attributes_2_0= ruleListAttribute ) ) ( (otherlv_3= ',' )? ( (lv_attributes_4_0= ruleListAttribute ) ) )* )? otherlv_5= ']'
                    {
                    otherlv_1=(Token)match(input,18,FOLLOW_12); 

                    				newLeafNode(otherlv_1, grammarAccess.getNodeStatementAccess().getLeftSquareBracketKeyword_1_0());
                    			
                    // InternalGraphvizDot.g:377:4: ( ( (lv_attributes_2_0= ruleListAttribute ) ) ( (otherlv_3= ',' )? ( (lv_attributes_4_0= ruleListAttribute ) ) )* )?
                    int alt9=2;
                    int LA9_0 = input.LA(1);

                    if ( ((LA9_0>=RULE_ID && LA9_0<=RULE_STRING)) ) {
                        alt9=1;
                    }
                    switch (alt9) {
                        case 1 :
                            // InternalGraphvizDot.g:378:5: ( (lv_attributes_2_0= ruleListAttribute ) ) ( (otherlv_3= ',' )? ( (lv_attributes_4_0= ruleListAttribute ) ) )*
                            {
                            // InternalGraphvizDot.g:378:5: ( (lv_attributes_2_0= ruleListAttribute ) )
                            // InternalGraphvizDot.g:379:6: (lv_attributes_2_0= ruleListAttribute )
                            {
                            // InternalGraphvizDot.g:379:6: (lv_attributes_2_0= ruleListAttribute )
                            // InternalGraphvizDot.g:380:7: lv_attributes_2_0= ruleListAttribute
                            {

                            							newCompositeNode(grammarAccess.getNodeStatementAccess().getAttributesListAttributeParserRuleCall_1_1_0_0());
                            						
                            pushFollow(FOLLOW_13);
                            lv_attributes_2_0=ruleListAttribute();

                            state._fsp--;


                            							if (current==null) {
                            								current = createModelElementForParent(grammarAccess.getNodeStatementRule());
                            							}
                            							add(
                            								current,
                            								"attributes",
                            								lv_attributes_2_0,
                            								"org.eclipse.elk.alg.graphviz.dot.GraphvizDot.ListAttribute");
                            							afterParserOrEnumRuleCall();
                            						

                            }


                            }

                            // InternalGraphvizDot.g:397:5: ( (otherlv_3= ',' )? ( (lv_attributes_4_0= ruleListAttribute ) ) )*
                            loop8:
                            do {
                                int alt8=2;
                                int LA8_0 = input.LA(1);

                                if ( ((LA8_0>=RULE_ID && LA8_0<=RULE_STRING)||LA8_0==19) ) {
                                    alt8=1;
                                }


                                switch (alt8) {
                            	case 1 :
                            	    // InternalGraphvizDot.g:398:6: (otherlv_3= ',' )? ( (lv_attributes_4_0= ruleListAttribute ) )
                            	    {
                            	    // InternalGraphvizDot.g:398:6: (otherlv_3= ',' )?
                            	    int alt7=2;
                            	    int LA7_0 = input.LA(1);

                            	    if ( (LA7_0==19) ) {
                            	        alt7=1;
                            	    }
                            	    switch (alt7) {
                            	        case 1 :
                            	            // InternalGraphvizDot.g:399:7: otherlv_3= ','
                            	            {
                            	            otherlv_3=(Token)match(input,19,FOLLOW_10); 

                            	            							newLeafNode(otherlv_3, grammarAccess.getNodeStatementAccess().getCommaKeyword_1_1_1_0());
                            	            						

                            	            }
                            	            break;

                            	    }

                            	    // InternalGraphvizDot.g:404:6: ( (lv_attributes_4_0= ruleListAttribute ) )
                            	    // InternalGraphvizDot.g:405:7: (lv_attributes_4_0= ruleListAttribute )
                            	    {
                            	    // InternalGraphvizDot.g:405:7: (lv_attributes_4_0= ruleListAttribute )
                            	    // InternalGraphvizDot.g:406:8: lv_attributes_4_0= ruleListAttribute
                            	    {

                            	    								newCompositeNode(grammarAccess.getNodeStatementAccess().getAttributesListAttributeParserRuleCall_1_1_1_1_0());
                            	    							
                            	    pushFollow(FOLLOW_13);
                            	    lv_attributes_4_0=ruleListAttribute();

                            	    state._fsp--;


                            	    								if (current==null) {
                            	    									current = createModelElementForParent(grammarAccess.getNodeStatementRule());
                            	    								}
                            	    								add(
                            	    									current,
                            	    									"attributes",
                            	    									lv_attributes_4_0,
                            	    									"org.eclipse.elk.alg.graphviz.dot.GraphvizDot.ListAttribute");
                            	    								afterParserOrEnumRuleCall();
                            	    							

                            	    }


                            	    }


                            	    }
                            	    break;

                            	default :
                            	    break loop8;
                                }
                            } while (true);


                            }
                            break;

                    }

                    otherlv_5=(Token)match(input,20,FOLLOW_2); 

                    				newLeafNode(otherlv_5, grammarAccess.getNodeStatementAccess().getRightSquareBracketKeyword_1_2());
                    			

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
    // $ANTLR end "ruleNodeStatement"


    // $ANTLR start "entryRuleNode"
    // InternalGraphvizDot.g:434:1: entryRuleNode returns [EObject current=null] : iv_ruleNode= ruleNode EOF ;
    public final EObject entryRuleNode() throws RecognitionException {
        EObject current = null;

        EObject iv_ruleNode = null;


        try {
            // InternalGraphvizDot.g:434:45: (iv_ruleNode= ruleNode EOF )
            // InternalGraphvizDot.g:435:2: iv_ruleNode= ruleNode EOF
            {
             newCompositeNode(grammarAccess.getNodeRule()); 
            pushFollow(FOLLOW_1);
            iv_ruleNode=ruleNode();

            state._fsp--;

             current =iv_ruleNode; 
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
    // $ANTLR end "entryRuleNode"


    // $ANTLR start "ruleNode"
    // InternalGraphvizDot.g:441:1: ruleNode returns [EObject current=null] : ( ( (lv_name_0_0= ruleDotID ) ) (otherlv_1= ':' ( (lv_port_2_0= rulePort ) ) )? ) ;
    public final EObject ruleNode() throws RecognitionException {
        EObject current = null;

        Token otherlv_1=null;
        AntlrDatatypeRuleToken lv_name_0_0 = null;

        EObject lv_port_2_0 = null;



        	enterRule();

        try {
            // InternalGraphvizDot.g:447:2: ( ( ( (lv_name_0_0= ruleDotID ) ) (otherlv_1= ':' ( (lv_port_2_0= rulePort ) ) )? ) )
            // InternalGraphvizDot.g:448:2: ( ( (lv_name_0_0= ruleDotID ) ) (otherlv_1= ':' ( (lv_port_2_0= rulePort ) ) )? )
            {
            // InternalGraphvizDot.g:448:2: ( ( (lv_name_0_0= ruleDotID ) ) (otherlv_1= ':' ( (lv_port_2_0= rulePort ) ) )? )
            // InternalGraphvizDot.g:449:3: ( (lv_name_0_0= ruleDotID ) ) (otherlv_1= ':' ( (lv_port_2_0= rulePort ) ) )?
            {
            // InternalGraphvizDot.g:449:3: ( (lv_name_0_0= ruleDotID ) )
            // InternalGraphvizDot.g:450:4: (lv_name_0_0= ruleDotID )
            {
            // InternalGraphvizDot.g:450:4: (lv_name_0_0= ruleDotID )
            // InternalGraphvizDot.g:451:5: lv_name_0_0= ruleDotID
            {

            					newCompositeNode(grammarAccess.getNodeAccess().getNameDotIDParserRuleCall_0_0());
            				
            pushFollow(FOLLOW_14);
            lv_name_0_0=ruleDotID();

            state._fsp--;


            					if (current==null) {
            						current = createModelElementForParent(grammarAccess.getNodeRule());
            					}
            					set(
            						current,
            						"name",
            						lv_name_0_0,
            						"org.eclipse.elk.alg.graphviz.dot.GraphvizDot.DotID");
            					afterParserOrEnumRuleCall();
            				

            }


            }

            // InternalGraphvizDot.g:468:3: (otherlv_1= ':' ( (lv_port_2_0= rulePort ) ) )?
            int alt11=2;
            int LA11_0 = input.LA(1);

            if ( (LA11_0==21) ) {
                alt11=1;
            }
            switch (alt11) {
                case 1 :
                    // InternalGraphvizDot.g:469:4: otherlv_1= ':' ( (lv_port_2_0= rulePort ) )
                    {
                    otherlv_1=(Token)match(input,21,FOLLOW_10); 

                    				newLeafNode(otherlv_1, grammarAccess.getNodeAccess().getColonKeyword_1_0());
                    			
                    // InternalGraphvizDot.g:473:4: ( (lv_port_2_0= rulePort ) )
                    // InternalGraphvizDot.g:474:5: (lv_port_2_0= rulePort )
                    {
                    // InternalGraphvizDot.g:474:5: (lv_port_2_0= rulePort )
                    // InternalGraphvizDot.g:475:6: lv_port_2_0= rulePort
                    {

                    						newCompositeNode(grammarAccess.getNodeAccess().getPortPortParserRuleCall_1_1_0());
                    					
                    pushFollow(FOLLOW_2);
                    lv_port_2_0=rulePort();

                    state._fsp--;


                    						if (current==null) {
                    							current = createModelElementForParent(grammarAccess.getNodeRule());
                    						}
                    						set(
                    							current,
                    							"port",
                    							lv_port_2_0,
                    							"org.eclipse.elk.alg.graphviz.dot.GraphvizDot.Port");
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
    // $ANTLR end "ruleNode"


    // $ANTLR start "entryRuleEdgeStatement"
    // InternalGraphvizDot.g:497:1: entryRuleEdgeStatement returns [EObject current=null] : iv_ruleEdgeStatement= ruleEdgeStatement EOF ;
    public final EObject entryRuleEdgeStatement() throws RecognitionException {
        EObject current = null;

        EObject iv_ruleEdgeStatement = null;


        try {
            // InternalGraphvizDot.g:497:54: (iv_ruleEdgeStatement= ruleEdgeStatement EOF )
            // InternalGraphvizDot.g:498:2: iv_ruleEdgeStatement= ruleEdgeStatement EOF
            {
             newCompositeNode(grammarAccess.getEdgeStatementRule()); 
            pushFollow(FOLLOW_1);
            iv_ruleEdgeStatement=ruleEdgeStatement();

            state._fsp--;

             current =iv_ruleEdgeStatement; 
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
    // $ANTLR end "entryRuleEdgeStatement"


    // $ANTLR start "ruleEdgeStatement"
    // InternalGraphvizDot.g:504:1: ruleEdgeStatement returns [EObject current=null] : ( ( (lv_sourceNode_0_0= ruleNode ) ) ( (lv_edgeTargets_1_0= ruleEdgeTarget ) )+ (otherlv_2= '[' ( ( (lv_attributes_3_0= ruleListAttribute ) ) ( (otherlv_4= ',' )? ( (lv_attributes_5_0= ruleListAttribute ) ) )* )? otherlv_6= ']' )? ) ;
    public final EObject ruleEdgeStatement() throws RecognitionException {
        EObject current = null;

        Token otherlv_2=null;
        Token otherlv_4=null;
        Token otherlv_6=null;
        EObject lv_sourceNode_0_0 = null;

        EObject lv_edgeTargets_1_0 = null;

        EObject lv_attributes_3_0 = null;

        EObject lv_attributes_5_0 = null;



        	enterRule();

        try {
            // InternalGraphvizDot.g:510:2: ( ( ( (lv_sourceNode_0_0= ruleNode ) ) ( (lv_edgeTargets_1_0= ruleEdgeTarget ) )+ (otherlv_2= '[' ( ( (lv_attributes_3_0= ruleListAttribute ) ) ( (otherlv_4= ',' )? ( (lv_attributes_5_0= ruleListAttribute ) ) )* )? otherlv_6= ']' )? ) )
            // InternalGraphvizDot.g:511:2: ( ( (lv_sourceNode_0_0= ruleNode ) ) ( (lv_edgeTargets_1_0= ruleEdgeTarget ) )+ (otherlv_2= '[' ( ( (lv_attributes_3_0= ruleListAttribute ) ) ( (otherlv_4= ',' )? ( (lv_attributes_5_0= ruleListAttribute ) ) )* )? otherlv_6= ']' )? )
            {
            // InternalGraphvizDot.g:511:2: ( ( (lv_sourceNode_0_0= ruleNode ) ) ( (lv_edgeTargets_1_0= ruleEdgeTarget ) )+ (otherlv_2= '[' ( ( (lv_attributes_3_0= ruleListAttribute ) ) ( (otherlv_4= ',' )? ( (lv_attributes_5_0= ruleListAttribute ) ) )* )? otherlv_6= ']' )? )
            // InternalGraphvizDot.g:512:3: ( (lv_sourceNode_0_0= ruleNode ) ) ( (lv_edgeTargets_1_0= ruleEdgeTarget ) )+ (otherlv_2= '[' ( ( (lv_attributes_3_0= ruleListAttribute ) ) ( (otherlv_4= ',' )? ( (lv_attributes_5_0= ruleListAttribute ) ) )* )? otherlv_6= ']' )?
            {
            // InternalGraphvizDot.g:512:3: ( (lv_sourceNode_0_0= ruleNode ) )
            // InternalGraphvizDot.g:513:4: (lv_sourceNode_0_0= ruleNode )
            {
            // InternalGraphvizDot.g:513:4: (lv_sourceNode_0_0= ruleNode )
            // InternalGraphvizDot.g:514:5: lv_sourceNode_0_0= ruleNode
            {

            					newCompositeNode(grammarAccess.getEdgeStatementAccess().getSourceNodeNodeParserRuleCall_0_0());
            				
            pushFollow(FOLLOW_15);
            lv_sourceNode_0_0=ruleNode();

            state._fsp--;


            					if (current==null) {
            						current = createModelElementForParent(grammarAccess.getEdgeStatementRule());
            					}
            					set(
            						current,
            						"sourceNode",
            						lv_sourceNode_0_0,
            						"org.eclipse.elk.alg.graphviz.dot.GraphvizDot.Node");
            					afterParserOrEnumRuleCall();
            				

            }


            }

            // InternalGraphvizDot.g:531:3: ( (lv_edgeTargets_1_0= ruleEdgeTarget ) )+
            int cnt12=0;
            loop12:
            do {
                int alt12=2;
                int LA12_0 = input.LA(1);

                if ( ((LA12_0>=23 && LA12_0<=24)) ) {
                    alt12=1;
                }


                switch (alt12) {
            	case 1 :
            	    // InternalGraphvizDot.g:532:4: (lv_edgeTargets_1_0= ruleEdgeTarget )
            	    {
            	    // InternalGraphvizDot.g:532:4: (lv_edgeTargets_1_0= ruleEdgeTarget )
            	    // InternalGraphvizDot.g:533:5: lv_edgeTargets_1_0= ruleEdgeTarget
            	    {

            	    					newCompositeNode(grammarAccess.getEdgeStatementAccess().getEdgeTargetsEdgeTargetParserRuleCall_1_0());
            	    				
            	    pushFollow(FOLLOW_16);
            	    lv_edgeTargets_1_0=ruleEdgeTarget();

            	    state._fsp--;


            	    					if (current==null) {
            	    						current = createModelElementForParent(grammarAccess.getEdgeStatementRule());
            	    					}
            	    					add(
            	    						current,
            	    						"edgeTargets",
            	    						lv_edgeTargets_1_0,
            	    						"org.eclipse.elk.alg.graphviz.dot.GraphvizDot.EdgeTarget");
            	    					afterParserOrEnumRuleCall();
            	    				

            	    }


            	    }
            	    break;

            	default :
            	    if ( cnt12 >= 1 ) break loop12;
                        EarlyExitException eee =
                            new EarlyExitException(12, input);
                        throw eee;
                }
                cnt12++;
            } while (true);

            // InternalGraphvizDot.g:550:3: (otherlv_2= '[' ( ( (lv_attributes_3_0= ruleListAttribute ) ) ( (otherlv_4= ',' )? ( (lv_attributes_5_0= ruleListAttribute ) ) )* )? otherlv_6= ']' )?
            int alt16=2;
            int LA16_0 = input.LA(1);

            if ( (LA16_0==18) ) {
                alt16=1;
            }
            switch (alt16) {
                case 1 :
                    // InternalGraphvizDot.g:551:4: otherlv_2= '[' ( ( (lv_attributes_3_0= ruleListAttribute ) ) ( (otherlv_4= ',' )? ( (lv_attributes_5_0= ruleListAttribute ) ) )* )? otherlv_6= ']'
                    {
                    otherlv_2=(Token)match(input,18,FOLLOW_12); 

                    				newLeafNode(otherlv_2, grammarAccess.getEdgeStatementAccess().getLeftSquareBracketKeyword_2_0());
                    			
                    // InternalGraphvizDot.g:555:4: ( ( (lv_attributes_3_0= ruleListAttribute ) ) ( (otherlv_4= ',' )? ( (lv_attributes_5_0= ruleListAttribute ) ) )* )?
                    int alt15=2;
                    int LA15_0 = input.LA(1);

                    if ( ((LA15_0>=RULE_ID && LA15_0<=RULE_STRING)) ) {
                        alt15=1;
                    }
                    switch (alt15) {
                        case 1 :
                            // InternalGraphvizDot.g:556:5: ( (lv_attributes_3_0= ruleListAttribute ) ) ( (otherlv_4= ',' )? ( (lv_attributes_5_0= ruleListAttribute ) ) )*
                            {
                            // InternalGraphvizDot.g:556:5: ( (lv_attributes_3_0= ruleListAttribute ) )
                            // InternalGraphvizDot.g:557:6: (lv_attributes_3_0= ruleListAttribute )
                            {
                            // InternalGraphvizDot.g:557:6: (lv_attributes_3_0= ruleListAttribute )
                            // InternalGraphvizDot.g:558:7: lv_attributes_3_0= ruleListAttribute
                            {

                            							newCompositeNode(grammarAccess.getEdgeStatementAccess().getAttributesListAttributeParserRuleCall_2_1_0_0());
                            						
                            pushFollow(FOLLOW_13);
                            lv_attributes_3_0=ruleListAttribute();

                            state._fsp--;


                            							if (current==null) {
                            								current = createModelElementForParent(grammarAccess.getEdgeStatementRule());
                            							}
                            							add(
                            								current,
                            								"attributes",
                            								lv_attributes_3_0,
                            								"org.eclipse.elk.alg.graphviz.dot.GraphvizDot.ListAttribute");
                            							afterParserOrEnumRuleCall();
                            						

                            }


                            }

                            // InternalGraphvizDot.g:575:5: ( (otherlv_4= ',' )? ( (lv_attributes_5_0= ruleListAttribute ) ) )*
                            loop14:
                            do {
                                int alt14=2;
                                int LA14_0 = input.LA(1);

                                if ( ((LA14_0>=RULE_ID && LA14_0<=RULE_STRING)||LA14_0==19) ) {
                                    alt14=1;
                                }


                                switch (alt14) {
                            	case 1 :
                            	    // InternalGraphvizDot.g:576:6: (otherlv_4= ',' )? ( (lv_attributes_5_0= ruleListAttribute ) )
                            	    {
                            	    // InternalGraphvizDot.g:576:6: (otherlv_4= ',' )?
                            	    int alt13=2;
                            	    int LA13_0 = input.LA(1);

                            	    if ( (LA13_0==19) ) {
                            	        alt13=1;
                            	    }
                            	    switch (alt13) {
                            	        case 1 :
                            	            // InternalGraphvizDot.g:577:7: otherlv_4= ','
                            	            {
                            	            otherlv_4=(Token)match(input,19,FOLLOW_10); 

                            	            							newLeafNode(otherlv_4, grammarAccess.getEdgeStatementAccess().getCommaKeyword_2_1_1_0());
                            	            						

                            	            }
                            	            break;

                            	    }

                            	    // InternalGraphvizDot.g:582:6: ( (lv_attributes_5_0= ruleListAttribute ) )
                            	    // InternalGraphvizDot.g:583:7: (lv_attributes_5_0= ruleListAttribute )
                            	    {
                            	    // InternalGraphvizDot.g:583:7: (lv_attributes_5_0= ruleListAttribute )
                            	    // InternalGraphvizDot.g:584:8: lv_attributes_5_0= ruleListAttribute
                            	    {

                            	    								newCompositeNode(grammarAccess.getEdgeStatementAccess().getAttributesListAttributeParserRuleCall_2_1_1_1_0());
                            	    							
                            	    pushFollow(FOLLOW_13);
                            	    lv_attributes_5_0=ruleListAttribute();

                            	    state._fsp--;


                            	    								if (current==null) {
                            	    									current = createModelElementForParent(grammarAccess.getEdgeStatementRule());
                            	    								}
                            	    								add(
                            	    									current,
                            	    									"attributes",
                            	    									lv_attributes_5_0,
                            	    									"org.eclipse.elk.alg.graphviz.dot.GraphvizDot.ListAttribute");
                            	    								afterParserOrEnumRuleCall();
                            	    							

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

                    otherlv_6=(Token)match(input,20,FOLLOW_2); 

                    				newLeafNode(otherlv_6, grammarAccess.getEdgeStatementAccess().getRightSquareBracketKeyword_2_2());
                    			

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
    // $ANTLR end "ruleEdgeStatement"


    // $ANTLR start "entryRuleEdgeTarget"
    // InternalGraphvizDot.g:612:1: entryRuleEdgeTarget returns [EObject current=null] : iv_ruleEdgeTarget= ruleEdgeTarget EOF ;
    public final EObject entryRuleEdgeTarget() throws RecognitionException {
        EObject current = null;

        EObject iv_ruleEdgeTarget = null;


        try {
            // InternalGraphvizDot.g:612:51: (iv_ruleEdgeTarget= ruleEdgeTarget EOF )
            // InternalGraphvizDot.g:613:2: iv_ruleEdgeTarget= ruleEdgeTarget EOF
            {
             newCompositeNode(grammarAccess.getEdgeTargetRule()); 
            pushFollow(FOLLOW_1);
            iv_ruleEdgeTarget=ruleEdgeTarget();

            state._fsp--;

             current =iv_ruleEdgeTarget; 
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
    // $ANTLR end "entryRuleEdgeTarget"


    // $ANTLR start "ruleEdgeTarget"
    // InternalGraphvizDot.g:619:1: ruleEdgeTarget returns [EObject current=null] : ( ( (lv_operator_0_0= ruleEdgeOperator ) ) ( ( (lv_targetSubgraph_1_0= ruleSubgraph ) ) | ( (lv_targetnode_2_0= ruleNode ) ) ) ) ;
    public final EObject ruleEdgeTarget() throws RecognitionException {
        EObject current = null;

        Enumerator lv_operator_0_0 = null;

        EObject lv_targetSubgraph_1_0 = null;

        EObject lv_targetnode_2_0 = null;



        	enterRule();

        try {
            // InternalGraphvizDot.g:625:2: ( ( ( (lv_operator_0_0= ruleEdgeOperator ) ) ( ( (lv_targetSubgraph_1_0= ruleSubgraph ) ) | ( (lv_targetnode_2_0= ruleNode ) ) ) ) )
            // InternalGraphvizDot.g:626:2: ( ( (lv_operator_0_0= ruleEdgeOperator ) ) ( ( (lv_targetSubgraph_1_0= ruleSubgraph ) ) | ( (lv_targetnode_2_0= ruleNode ) ) ) )
            {
            // InternalGraphvizDot.g:626:2: ( ( (lv_operator_0_0= ruleEdgeOperator ) ) ( ( (lv_targetSubgraph_1_0= ruleSubgraph ) ) | ( (lv_targetnode_2_0= ruleNode ) ) ) )
            // InternalGraphvizDot.g:627:3: ( (lv_operator_0_0= ruleEdgeOperator ) ) ( ( (lv_targetSubgraph_1_0= ruleSubgraph ) ) | ( (lv_targetnode_2_0= ruleNode ) ) )
            {
            // InternalGraphvizDot.g:627:3: ( (lv_operator_0_0= ruleEdgeOperator ) )
            // InternalGraphvizDot.g:628:4: (lv_operator_0_0= ruleEdgeOperator )
            {
            // InternalGraphvizDot.g:628:4: (lv_operator_0_0= ruleEdgeOperator )
            // InternalGraphvizDot.g:629:5: lv_operator_0_0= ruleEdgeOperator
            {

            					newCompositeNode(grammarAccess.getEdgeTargetAccess().getOperatorEdgeOperatorEnumRuleCall_0_0());
            				
            pushFollow(FOLLOW_17);
            lv_operator_0_0=ruleEdgeOperator();

            state._fsp--;


            					if (current==null) {
            						current = createModelElementForParent(grammarAccess.getEdgeTargetRule());
            					}
            					set(
            						current,
            						"operator",
            						lv_operator_0_0,
            						"org.eclipse.elk.alg.graphviz.dot.GraphvizDot.EdgeOperator");
            					afterParserOrEnumRuleCall();
            				

            }


            }

            // InternalGraphvizDot.g:646:3: ( ( (lv_targetSubgraph_1_0= ruleSubgraph ) ) | ( (lv_targetnode_2_0= ruleNode ) ) )
            int alt17=2;
            int LA17_0 = input.LA(1);

            if ( (LA17_0==14||LA17_0==22) ) {
                alt17=1;
            }
            else if ( ((LA17_0>=RULE_ID && LA17_0<=RULE_STRING)) ) {
                alt17=2;
            }
            else {
                NoViableAltException nvae =
                    new NoViableAltException("", 17, 0, input);

                throw nvae;
            }
            switch (alt17) {
                case 1 :
                    // InternalGraphvizDot.g:647:4: ( (lv_targetSubgraph_1_0= ruleSubgraph ) )
                    {
                    // InternalGraphvizDot.g:647:4: ( (lv_targetSubgraph_1_0= ruleSubgraph ) )
                    // InternalGraphvizDot.g:648:5: (lv_targetSubgraph_1_0= ruleSubgraph )
                    {
                    // InternalGraphvizDot.g:648:5: (lv_targetSubgraph_1_0= ruleSubgraph )
                    // InternalGraphvizDot.g:649:6: lv_targetSubgraph_1_0= ruleSubgraph
                    {

                    						newCompositeNode(grammarAccess.getEdgeTargetAccess().getTargetSubgraphSubgraphParserRuleCall_1_0_0());
                    					
                    pushFollow(FOLLOW_2);
                    lv_targetSubgraph_1_0=ruleSubgraph();

                    state._fsp--;


                    						if (current==null) {
                    							current = createModelElementForParent(grammarAccess.getEdgeTargetRule());
                    						}
                    						set(
                    							current,
                    							"targetSubgraph",
                    							lv_targetSubgraph_1_0,
                    							"org.eclipse.elk.alg.graphviz.dot.GraphvizDot.Subgraph");
                    						afterParserOrEnumRuleCall();
                    					

                    }


                    }


                    }
                    break;
                case 2 :
                    // InternalGraphvizDot.g:667:4: ( (lv_targetnode_2_0= ruleNode ) )
                    {
                    // InternalGraphvizDot.g:667:4: ( (lv_targetnode_2_0= ruleNode ) )
                    // InternalGraphvizDot.g:668:5: (lv_targetnode_2_0= ruleNode )
                    {
                    // InternalGraphvizDot.g:668:5: (lv_targetnode_2_0= ruleNode )
                    // InternalGraphvizDot.g:669:6: lv_targetnode_2_0= ruleNode
                    {

                    						newCompositeNode(grammarAccess.getEdgeTargetAccess().getTargetnodeNodeParserRuleCall_1_1_0());
                    					
                    pushFollow(FOLLOW_2);
                    lv_targetnode_2_0=ruleNode();

                    state._fsp--;


                    						if (current==null) {
                    							current = createModelElementForParent(grammarAccess.getEdgeTargetRule());
                    						}
                    						set(
                    							current,
                    							"targetnode",
                    							lv_targetnode_2_0,
                    							"org.eclipse.elk.alg.graphviz.dot.GraphvizDot.Node");
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
    // $ANTLR end "ruleEdgeTarget"


    // $ANTLR start "entryRuleAttributeStatement"
    // InternalGraphvizDot.g:691:1: entryRuleAttributeStatement returns [EObject current=null] : iv_ruleAttributeStatement= ruleAttributeStatement EOF ;
    public final EObject entryRuleAttributeStatement() throws RecognitionException {
        EObject current = null;

        EObject iv_ruleAttributeStatement = null;


        try {
            // InternalGraphvizDot.g:691:59: (iv_ruleAttributeStatement= ruleAttributeStatement EOF )
            // InternalGraphvizDot.g:692:2: iv_ruleAttributeStatement= ruleAttributeStatement EOF
            {
             newCompositeNode(grammarAccess.getAttributeStatementRule()); 
            pushFollow(FOLLOW_1);
            iv_ruleAttributeStatement=ruleAttributeStatement();

            state._fsp--;

             current =iv_ruleAttributeStatement; 
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
    // $ANTLR end "entryRuleAttributeStatement"


    // $ANTLR start "ruleAttributeStatement"
    // InternalGraphvizDot.g:698:1: ruleAttributeStatement returns [EObject current=null] : ( ( (lv_type_0_0= ruleAttributeType ) ) otherlv_1= '[' ( ( (lv_attributes_2_0= ruleListAttribute ) ) ( (otherlv_3= ',' )? ( (lv_attributes_4_0= ruleListAttribute ) ) )* )? otherlv_5= ']' ) ;
    public final EObject ruleAttributeStatement() throws RecognitionException {
        EObject current = null;

        Token otherlv_1=null;
        Token otherlv_3=null;
        Token otherlv_5=null;
        Enumerator lv_type_0_0 = null;

        EObject lv_attributes_2_0 = null;

        EObject lv_attributes_4_0 = null;



        	enterRule();

        try {
            // InternalGraphvizDot.g:704:2: ( ( ( (lv_type_0_0= ruleAttributeType ) ) otherlv_1= '[' ( ( (lv_attributes_2_0= ruleListAttribute ) ) ( (otherlv_3= ',' )? ( (lv_attributes_4_0= ruleListAttribute ) ) )* )? otherlv_5= ']' ) )
            // InternalGraphvizDot.g:705:2: ( ( (lv_type_0_0= ruleAttributeType ) ) otherlv_1= '[' ( ( (lv_attributes_2_0= ruleListAttribute ) ) ( (otherlv_3= ',' )? ( (lv_attributes_4_0= ruleListAttribute ) ) )* )? otherlv_5= ']' )
            {
            // InternalGraphvizDot.g:705:2: ( ( (lv_type_0_0= ruleAttributeType ) ) otherlv_1= '[' ( ( (lv_attributes_2_0= ruleListAttribute ) ) ( (otherlv_3= ',' )? ( (lv_attributes_4_0= ruleListAttribute ) ) )* )? otherlv_5= ']' )
            // InternalGraphvizDot.g:706:3: ( (lv_type_0_0= ruleAttributeType ) ) otherlv_1= '[' ( ( (lv_attributes_2_0= ruleListAttribute ) ) ( (otherlv_3= ',' )? ( (lv_attributes_4_0= ruleListAttribute ) ) )* )? otherlv_5= ']'
            {
            // InternalGraphvizDot.g:706:3: ( (lv_type_0_0= ruleAttributeType ) )
            // InternalGraphvizDot.g:707:4: (lv_type_0_0= ruleAttributeType )
            {
            // InternalGraphvizDot.g:707:4: (lv_type_0_0= ruleAttributeType )
            // InternalGraphvizDot.g:708:5: lv_type_0_0= ruleAttributeType
            {

            					newCompositeNode(grammarAccess.getAttributeStatementAccess().getTypeAttributeTypeEnumRuleCall_0_0());
            				
            pushFollow(FOLLOW_18);
            lv_type_0_0=ruleAttributeType();

            state._fsp--;


            					if (current==null) {
            						current = createModelElementForParent(grammarAccess.getAttributeStatementRule());
            					}
            					set(
            						current,
            						"type",
            						lv_type_0_0,
            						"org.eclipse.elk.alg.graphviz.dot.GraphvizDot.AttributeType");
            					afterParserOrEnumRuleCall();
            				

            }


            }

            otherlv_1=(Token)match(input,18,FOLLOW_12); 

            			newLeafNode(otherlv_1, grammarAccess.getAttributeStatementAccess().getLeftSquareBracketKeyword_1());
            		
            // InternalGraphvizDot.g:729:3: ( ( (lv_attributes_2_0= ruleListAttribute ) ) ( (otherlv_3= ',' )? ( (lv_attributes_4_0= ruleListAttribute ) ) )* )?
            int alt20=2;
            int LA20_0 = input.LA(1);

            if ( ((LA20_0>=RULE_ID && LA20_0<=RULE_STRING)) ) {
                alt20=1;
            }
            switch (alt20) {
                case 1 :
                    // InternalGraphvizDot.g:730:4: ( (lv_attributes_2_0= ruleListAttribute ) ) ( (otherlv_3= ',' )? ( (lv_attributes_4_0= ruleListAttribute ) ) )*
                    {
                    // InternalGraphvizDot.g:730:4: ( (lv_attributes_2_0= ruleListAttribute ) )
                    // InternalGraphvizDot.g:731:5: (lv_attributes_2_0= ruleListAttribute )
                    {
                    // InternalGraphvizDot.g:731:5: (lv_attributes_2_0= ruleListAttribute )
                    // InternalGraphvizDot.g:732:6: lv_attributes_2_0= ruleListAttribute
                    {

                    						newCompositeNode(grammarAccess.getAttributeStatementAccess().getAttributesListAttributeParserRuleCall_2_0_0());
                    					
                    pushFollow(FOLLOW_13);
                    lv_attributes_2_0=ruleListAttribute();

                    state._fsp--;


                    						if (current==null) {
                    							current = createModelElementForParent(grammarAccess.getAttributeStatementRule());
                    						}
                    						add(
                    							current,
                    							"attributes",
                    							lv_attributes_2_0,
                    							"org.eclipse.elk.alg.graphviz.dot.GraphvizDot.ListAttribute");
                    						afterParserOrEnumRuleCall();
                    					

                    }


                    }

                    // InternalGraphvizDot.g:749:4: ( (otherlv_3= ',' )? ( (lv_attributes_4_0= ruleListAttribute ) ) )*
                    loop19:
                    do {
                        int alt19=2;
                        int LA19_0 = input.LA(1);

                        if ( ((LA19_0>=RULE_ID && LA19_0<=RULE_STRING)||LA19_0==19) ) {
                            alt19=1;
                        }


                        switch (alt19) {
                    	case 1 :
                    	    // InternalGraphvizDot.g:750:5: (otherlv_3= ',' )? ( (lv_attributes_4_0= ruleListAttribute ) )
                    	    {
                    	    // InternalGraphvizDot.g:750:5: (otherlv_3= ',' )?
                    	    int alt18=2;
                    	    int LA18_0 = input.LA(1);

                    	    if ( (LA18_0==19) ) {
                    	        alt18=1;
                    	    }
                    	    switch (alt18) {
                    	        case 1 :
                    	            // InternalGraphvizDot.g:751:6: otherlv_3= ','
                    	            {
                    	            otherlv_3=(Token)match(input,19,FOLLOW_10); 

                    	            						newLeafNode(otherlv_3, grammarAccess.getAttributeStatementAccess().getCommaKeyword_2_1_0());
                    	            					

                    	            }
                    	            break;

                    	    }

                    	    // InternalGraphvizDot.g:756:5: ( (lv_attributes_4_0= ruleListAttribute ) )
                    	    // InternalGraphvizDot.g:757:6: (lv_attributes_4_0= ruleListAttribute )
                    	    {
                    	    // InternalGraphvizDot.g:757:6: (lv_attributes_4_0= ruleListAttribute )
                    	    // InternalGraphvizDot.g:758:7: lv_attributes_4_0= ruleListAttribute
                    	    {

                    	    							newCompositeNode(grammarAccess.getAttributeStatementAccess().getAttributesListAttributeParserRuleCall_2_1_1_0());
                    	    						
                    	    pushFollow(FOLLOW_13);
                    	    lv_attributes_4_0=ruleListAttribute();

                    	    state._fsp--;


                    	    							if (current==null) {
                    	    								current = createModelElementForParent(grammarAccess.getAttributeStatementRule());
                    	    							}
                    	    							add(
                    	    								current,
                    	    								"attributes",
                    	    								lv_attributes_4_0,
                    	    								"org.eclipse.elk.alg.graphviz.dot.GraphvizDot.ListAttribute");
                    	    							afterParserOrEnumRuleCall();
                    	    						

                    	    }


                    	    }


                    	    }
                    	    break;

                    	default :
                    	    break loop19;
                        }
                    } while (true);


                    }
                    break;

            }

            otherlv_5=(Token)match(input,20,FOLLOW_2); 

            			newLeafNode(otherlv_5, grammarAccess.getAttributeStatementAccess().getRightSquareBracketKeyword_3());
            		

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
    // $ANTLR end "ruleAttributeStatement"


    // $ANTLR start "entryRuleSubgraph"
    // InternalGraphvizDot.g:785:1: entryRuleSubgraph returns [EObject current=null] : iv_ruleSubgraph= ruleSubgraph EOF ;
    public final EObject entryRuleSubgraph() throws RecognitionException {
        EObject current = null;

        EObject iv_ruleSubgraph = null;


        try {
            // InternalGraphvizDot.g:785:49: (iv_ruleSubgraph= ruleSubgraph EOF )
            // InternalGraphvizDot.g:786:2: iv_ruleSubgraph= ruleSubgraph EOF
            {
             newCompositeNode(grammarAccess.getSubgraphRule()); 
            pushFollow(FOLLOW_1);
            iv_ruleSubgraph=ruleSubgraph();

            state._fsp--;

             current =iv_ruleSubgraph; 
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
    // $ANTLR end "entryRuleSubgraph"


    // $ANTLR start "ruleSubgraph"
    // InternalGraphvizDot.g:792:1: ruleSubgraph returns [EObject current=null] : ( () (otherlv_1= 'subgraph' ( (lv_name_2_0= RULE_ID ) )? )? otherlv_3= '{' ( (lv_statements_4_0= ruleStatement ) )* otherlv_5= '}' ) ;
    public final EObject ruleSubgraph() throws RecognitionException {
        EObject current = null;

        Token otherlv_1=null;
        Token lv_name_2_0=null;
        Token otherlv_3=null;
        Token otherlv_5=null;
        EObject lv_statements_4_0 = null;



        	enterRule();

        try {
            // InternalGraphvizDot.g:798:2: ( ( () (otherlv_1= 'subgraph' ( (lv_name_2_0= RULE_ID ) )? )? otherlv_3= '{' ( (lv_statements_4_0= ruleStatement ) )* otherlv_5= '}' ) )
            // InternalGraphvizDot.g:799:2: ( () (otherlv_1= 'subgraph' ( (lv_name_2_0= RULE_ID ) )? )? otherlv_3= '{' ( (lv_statements_4_0= ruleStatement ) )* otherlv_5= '}' )
            {
            // InternalGraphvizDot.g:799:2: ( () (otherlv_1= 'subgraph' ( (lv_name_2_0= RULE_ID ) )? )? otherlv_3= '{' ( (lv_statements_4_0= ruleStatement ) )* otherlv_5= '}' )
            // InternalGraphvizDot.g:800:3: () (otherlv_1= 'subgraph' ( (lv_name_2_0= RULE_ID ) )? )? otherlv_3= '{' ( (lv_statements_4_0= ruleStatement ) )* otherlv_5= '}'
            {
            // InternalGraphvizDot.g:800:3: ()
            // InternalGraphvizDot.g:801:4: 
            {

            				current = forceCreateModelElement(
            					grammarAccess.getSubgraphAccess().getSubgraphAction_0(),
            					current);
            			

            }

            // InternalGraphvizDot.g:807:3: (otherlv_1= 'subgraph' ( (lv_name_2_0= RULE_ID ) )? )?
            int alt22=2;
            int LA22_0 = input.LA(1);

            if ( (LA22_0==22) ) {
                alt22=1;
            }
            switch (alt22) {
                case 1 :
                    // InternalGraphvizDot.g:808:4: otherlv_1= 'subgraph' ( (lv_name_2_0= RULE_ID ) )?
                    {
                    otherlv_1=(Token)match(input,22,FOLLOW_19); 

                    				newLeafNode(otherlv_1, grammarAccess.getSubgraphAccess().getSubgraphKeyword_1_0());
                    			
                    // InternalGraphvizDot.g:812:4: ( (lv_name_2_0= RULE_ID ) )?
                    int alt21=2;
                    int LA21_0 = input.LA(1);

                    if ( (LA21_0==RULE_ID) ) {
                        alt21=1;
                    }
                    switch (alt21) {
                        case 1 :
                            // InternalGraphvizDot.g:813:5: (lv_name_2_0= RULE_ID )
                            {
                            // InternalGraphvizDot.g:813:5: (lv_name_2_0= RULE_ID )
                            // InternalGraphvizDot.g:814:6: lv_name_2_0= RULE_ID
                            {
                            lv_name_2_0=(Token)match(input,RULE_ID,FOLLOW_6); 

                            						newLeafNode(lv_name_2_0, grammarAccess.getSubgraphAccess().getNameIDTerminalRuleCall_1_1_0());
                            					

                            						if (current==null) {
                            							current = createModelElement(grammarAccess.getSubgraphRule());
                            						}
                            						setWithLastConsumed(
                            							current,
                            							"name",
                            							lv_name_2_0,
                            							"org.eclipse.elk.alg.graphviz.dot.GraphvizDot.ID");
                            					

                            }


                            }
                            break;

                    }


                    }
                    break;

            }

            otherlv_3=(Token)match(input,14,FOLLOW_7); 

            			newLeafNode(otherlv_3, grammarAccess.getSubgraphAccess().getLeftCurlyBracketKeyword_2());
            		
            // InternalGraphvizDot.g:835:3: ( (lv_statements_4_0= ruleStatement ) )*
            loop23:
            do {
                int alt23=2;
                int LA23_0 = input.LA(1);

                if ( ((LA23_0>=RULE_ID && LA23_0<=RULE_STRING)||LA23_0==14||LA23_0==22||LA23_0==25||(LA23_0>=27 && LA23_0<=28)) ) {
                    alt23=1;
                }


                switch (alt23) {
            	case 1 :
            	    // InternalGraphvizDot.g:836:4: (lv_statements_4_0= ruleStatement )
            	    {
            	    // InternalGraphvizDot.g:836:4: (lv_statements_4_0= ruleStatement )
            	    // InternalGraphvizDot.g:837:5: lv_statements_4_0= ruleStatement
            	    {

            	    					newCompositeNode(grammarAccess.getSubgraphAccess().getStatementsStatementParserRuleCall_3_0());
            	    				
            	    pushFollow(FOLLOW_7);
            	    lv_statements_4_0=ruleStatement();

            	    state._fsp--;


            	    					if (current==null) {
            	    						current = createModelElementForParent(grammarAccess.getSubgraphRule());
            	    					}
            	    					add(
            	    						current,
            	    						"statements",
            	    						lv_statements_4_0,
            	    						"org.eclipse.elk.alg.graphviz.dot.GraphvizDot.Statement");
            	    					afterParserOrEnumRuleCall();
            	    				

            	    }


            	    }
            	    break;

            	default :
            	    break loop23;
                }
            } while (true);

            otherlv_5=(Token)match(input,15,FOLLOW_2); 

            			newLeafNode(otherlv_5, grammarAccess.getSubgraphAccess().getRightCurlyBracketKeyword_4());
            		

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
    // $ANTLR end "ruleSubgraph"


    // $ANTLR start "entryRuleListAttribute"
    // InternalGraphvizDot.g:862:1: entryRuleListAttribute returns [EObject current=null] : iv_ruleListAttribute= ruleListAttribute EOF ;
    public final EObject entryRuleListAttribute() throws RecognitionException {
        EObject current = null;

        EObject iv_ruleListAttribute = null;


        try {
            // InternalGraphvizDot.g:862:54: (iv_ruleListAttribute= ruleListAttribute EOF )
            // InternalGraphvizDot.g:863:2: iv_ruleListAttribute= ruleListAttribute EOF
            {
             newCompositeNode(grammarAccess.getListAttributeRule()); 
            pushFollow(FOLLOW_1);
            iv_ruleListAttribute=ruleListAttribute();

            state._fsp--;

             current =iv_ruleListAttribute; 
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
    // $ANTLR end "entryRuleListAttribute"


    // $ANTLR start "ruleListAttribute"
    // InternalGraphvizDot.g:869:1: ruleListAttribute returns [EObject current=null] : ( ( (lv_name_0_0= ruleDotID ) ) (otherlv_1= '=' ( (lv_value_2_0= ruleDotID ) ) )? ) ;
    public final EObject ruleListAttribute() throws RecognitionException {
        EObject current = null;

        Token otherlv_1=null;
        AntlrDatatypeRuleToken lv_name_0_0 = null;

        AntlrDatatypeRuleToken lv_value_2_0 = null;



        	enterRule();

        try {
            // InternalGraphvizDot.g:875:2: ( ( ( (lv_name_0_0= ruleDotID ) ) (otherlv_1= '=' ( (lv_value_2_0= ruleDotID ) ) )? ) )
            // InternalGraphvizDot.g:876:2: ( ( (lv_name_0_0= ruleDotID ) ) (otherlv_1= '=' ( (lv_value_2_0= ruleDotID ) ) )? )
            {
            // InternalGraphvizDot.g:876:2: ( ( (lv_name_0_0= ruleDotID ) ) (otherlv_1= '=' ( (lv_value_2_0= ruleDotID ) ) )? )
            // InternalGraphvizDot.g:877:3: ( (lv_name_0_0= ruleDotID ) ) (otherlv_1= '=' ( (lv_value_2_0= ruleDotID ) ) )?
            {
            // InternalGraphvizDot.g:877:3: ( (lv_name_0_0= ruleDotID ) )
            // InternalGraphvizDot.g:878:4: (lv_name_0_0= ruleDotID )
            {
            // InternalGraphvizDot.g:878:4: (lv_name_0_0= ruleDotID )
            // InternalGraphvizDot.g:879:5: lv_name_0_0= ruleDotID
            {

            					newCompositeNode(grammarAccess.getListAttributeAccess().getNameDotIDParserRuleCall_0_0());
            				
            pushFollow(FOLLOW_20);
            lv_name_0_0=ruleDotID();

            state._fsp--;


            					if (current==null) {
            						current = createModelElementForParent(grammarAccess.getListAttributeRule());
            					}
            					set(
            						current,
            						"name",
            						lv_name_0_0,
            						"org.eclipse.elk.alg.graphviz.dot.GraphvizDot.DotID");
            					afterParserOrEnumRuleCall();
            				

            }


            }

            // InternalGraphvizDot.g:896:3: (otherlv_1= '=' ( (lv_value_2_0= ruleDotID ) ) )?
            int alt24=2;
            int LA24_0 = input.LA(1);

            if ( (LA24_0==17) ) {
                alt24=1;
            }
            switch (alt24) {
                case 1 :
                    // InternalGraphvizDot.g:897:4: otherlv_1= '=' ( (lv_value_2_0= ruleDotID ) )
                    {
                    otherlv_1=(Token)match(input,17,FOLLOW_10); 

                    				newLeafNode(otherlv_1, grammarAccess.getListAttributeAccess().getEqualsSignKeyword_1_0());
                    			
                    // InternalGraphvizDot.g:901:4: ( (lv_value_2_0= ruleDotID ) )
                    // InternalGraphvizDot.g:902:5: (lv_value_2_0= ruleDotID )
                    {
                    // InternalGraphvizDot.g:902:5: (lv_value_2_0= ruleDotID )
                    // InternalGraphvizDot.g:903:6: lv_value_2_0= ruleDotID
                    {

                    						newCompositeNode(grammarAccess.getListAttributeAccess().getValueDotIDParserRuleCall_1_1_0());
                    					
                    pushFollow(FOLLOW_2);
                    lv_value_2_0=ruleDotID();

                    state._fsp--;


                    						if (current==null) {
                    							current = createModelElementForParent(grammarAccess.getListAttributeRule());
                    						}
                    						set(
                    							current,
                    							"value",
                    							lv_value_2_0,
                    							"org.eclipse.elk.alg.graphviz.dot.GraphvizDot.DotID");
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
    // $ANTLR end "ruleListAttribute"


    // $ANTLR start "entryRulePort"
    // InternalGraphvizDot.g:925:1: entryRulePort returns [EObject current=null] : iv_rulePort= rulePort EOF ;
    public final EObject entryRulePort() throws RecognitionException {
        EObject current = null;

        EObject iv_rulePort = null;


        try {
            // InternalGraphvizDot.g:925:45: (iv_rulePort= rulePort EOF )
            // InternalGraphvizDot.g:926:2: iv_rulePort= rulePort EOF
            {
             newCompositeNode(grammarAccess.getPortRule()); 
            pushFollow(FOLLOW_1);
            iv_rulePort=rulePort();

            state._fsp--;

             current =iv_rulePort; 
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
    // $ANTLR end "entryRulePort"


    // $ANTLR start "rulePort"
    // InternalGraphvizDot.g:932:1: rulePort returns [EObject current=null] : ( ( (lv_name_0_0= ruleDotID ) ) (otherlv_1= ':' ( (lv_compass_pt_2_0= RULE_ID ) ) )? ) ;
    public final EObject rulePort() throws RecognitionException {
        EObject current = null;

        Token otherlv_1=null;
        Token lv_compass_pt_2_0=null;
        AntlrDatatypeRuleToken lv_name_0_0 = null;



        	enterRule();

        try {
            // InternalGraphvizDot.g:938:2: ( ( ( (lv_name_0_0= ruleDotID ) ) (otherlv_1= ':' ( (lv_compass_pt_2_0= RULE_ID ) ) )? ) )
            // InternalGraphvizDot.g:939:2: ( ( (lv_name_0_0= ruleDotID ) ) (otherlv_1= ':' ( (lv_compass_pt_2_0= RULE_ID ) ) )? )
            {
            // InternalGraphvizDot.g:939:2: ( ( (lv_name_0_0= ruleDotID ) ) (otherlv_1= ':' ( (lv_compass_pt_2_0= RULE_ID ) ) )? )
            // InternalGraphvizDot.g:940:3: ( (lv_name_0_0= ruleDotID ) ) (otherlv_1= ':' ( (lv_compass_pt_2_0= RULE_ID ) ) )?
            {
            // InternalGraphvizDot.g:940:3: ( (lv_name_0_0= ruleDotID ) )
            // InternalGraphvizDot.g:941:4: (lv_name_0_0= ruleDotID )
            {
            // InternalGraphvizDot.g:941:4: (lv_name_0_0= ruleDotID )
            // InternalGraphvizDot.g:942:5: lv_name_0_0= ruleDotID
            {

            					newCompositeNode(grammarAccess.getPortAccess().getNameDotIDParserRuleCall_0_0());
            				
            pushFollow(FOLLOW_14);
            lv_name_0_0=ruleDotID();

            state._fsp--;


            					if (current==null) {
            						current = createModelElementForParent(grammarAccess.getPortRule());
            					}
            					set(
            						current,
            						"name",
            						lv_name_0_0,
            						"org.eclipse.elk.alg.graphviz.dot.GraphvizDot.DotID");
            					afterParserOrEnumRuleCall();
            				

            }


            }

            // InternalGraphvizDot.g:959:3: (otherlv_1= ':' ( (lv_compass_pt_2_0= RULE_ID ) ) )?
            int alt25=2;
            int LA25_0 = input.LA(1);

            if ( (LA25_0==21) ) {
                alt25=1;
            }
            switch (alt25) {
                case 1 :
                    // InternalGraphvizDot.g:960:4: otherlv_1= ':' ( (lv_compass_pt_2_0= RULE_ID ) )
                    {
                    otherlv_1=(Token)match(input,21,FOLLOW_21); 

                    				newLeafNode(otherlv_1, grammarAccess.getPortAccess().getColonKeyword_1_0());
                    			
                    // InternalGraphvizDot.g:964:4: ( (lv_compass_pt_2_0= RULE_ID ) )
                    // InternalGraphvizDot.g:965:5: (lv_compass_pt_2_0= RULE_ID )
                    {
                    // InternalGraphvizDot.g:965:5: (lv_compass_pt_2_0= RULE_ID )
                    // InternalGraphvizDot.g:966:6: lv_compass_pt_2_0= RULE_ID
                    {
                    lv_compass_pt_2_0=(Token)match(input,RULE_ID,FOLLOW_2); 

                    						newLeafNode(lv_compass_pt_2_0, grammarAccess.getPortAccess().getCompass_ptIDTerminalRuleCall_1_1_0());
                    					

                    						if (current==null) {
                    							current = createModelElement(grammarAccess.getPortRule());
                    						}
                    						setWithLastConsumed(
                    							current,
                    							"compass_pt",
                    							lv_compass_pt_2_0,
                    							"org.eclipse.elk.alg.graphviz.dot.GraphvizDot.ID");
                    					

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
    // $ANTLR end "rulePort"


    // $ANTLR start "entryRuleDotID"
    // InternalGraphvizDot.g:987:1: entryRuleDotID returns [String current=null] : iv_ruleDotID= ruleDotID EOF ;
    public final String entryRuleDotID() throws RecognitionException {
        String current = null;

        AntlrDatatypeRuleToken iv_ruleDotID = null;


        try {
            // InternalGraphvizDot.g:987:45: (iv_ruleDotID= ruleDotID EOF )
            // InternalGraphvizDot.g:988:2: iv_ruleDotID= ruleDotID EOF
            {
             newCompositeNode(grammarAccess.getDotIDRule()); 
            pushFollow(FOLLOW_1);
            iv_ruleDotID=ruleDotID();

            state._fsp--;

             current =iv_ruleDotID.getText(); 
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
    // $ANTLR end "entryRuleDotID"


    // $ANTLR start "ruleDotID"
    // InternalGraphvizDot.g:994:1: ruleDotID returns [AntlrDatatypeRuleToken current=new AntlrDatatypeRuleToken()] : (this_ID_0= RULE_ID | this_INT_1= RULE_INT | this_FLOAT_2= RULE_FLOAT | this_STRING_3= RULE_STRING ) ;
    public final AntlrDatatypeRuleToken ruleDotID() throws RecognitionException {
        AntlrDatatypeRuleToken current = new AntlrDatatypeRuleToken();

        Token this_ID_0=null;
        Token this_INT_1=null;
        Token this_FLOAT_2=null;
        Token this_STRING_3=null;


        	enterRule();

        try {
            // InternalGraphvizDot.g:1000:2: ( (this_ID_0= RULE_ID | this_INT_1= RULE_INT | this_FLOAT_2= RULE_FLOAT | this_STRING_3= RULE_STRING ) )
            // InternalGraphvizDot.g:1001:2: (this_ID_0= RULE_ID | this_INT_1= RULE_INT | this_FLOAT_2= RULE_FLOAT | this_STRING_3= RULE_STRING )
            {
            // InternalGraphvizDot.g:1001:2: (this_ID_0= RULE_ID | this_INT_1= RULE_INT | this_FLOAT_2= RULE_FLOAT | this_STRING_3= RULE_STRING )
            int alt26=4;
            switch ( input.LA(1) ) {
            case RULE_ID:
                {
                alt26=1;
                }
                break;
            case RULE_INT:
                {
                alt26=2;
                }
                break;
            case RULE_FLOAT:
                {
                alt26=3;
                }
                break;
            case RULE_STRING:
                {
                alt26=4;
                }
                break;
            default:
                NoViableAltException nvae =
                    new NoViableAltException("", 26, 0, input);

                throw nvae;
            }

            switch (alt26) {
                case 1 :
                    // InternalGraphvizDot.g:1002:3: this_ID_0= RULE_ID
                    {
                    this_ID_0=(Token)match(input,RULE_ID,FOLLOW_2); 

                    			current.merge(this_ID_0);
                    		

                    			newLeafNode(this_ID_0, grammarAccess.getDotIDAccess().getIDTerminalRuleCall_0());
                    		

                    }
                    break;
                case 2 :
                    // InternalGraphvizDot.g:1010:3: this_INT_1= RULE_INT
                    {
                    this_INT_1=(Token)match(input,RULE_INT,FOLLOW_2); 

                    			current.merge(this_INT_1);
                    		

                    			newLeafNode(this_INT_1, grammarAccess.getDotIDAccess().getINTTerminalRuleCall_1());
                    		

                    }
                    break;
                case 3 :
                    // InternalGraphvizDot.g:1018:3: this_FLOAT_2= RULE_FLOAT
                    {
                    this_FLOAT_2=(Token)match(input,RULE_FLOAT,FOLLOW_2); 

                    			current.merge(this_FLOAT_2);
                    		

                    			newLeafNode(this_FLOAT_2, grammarAccess.getDotIDAccess().getFLOATTerminalRuleCall_2());
                    		

                    }
                    break;
                case 4 :
                    // InternalGraphvizDot.g:1026:3: this_STRING_3= RULE_STRING
                    {
                    this_STRING_3=(Token)match(input,RULE_STRING,FOLLOW_2); 

                    			current.merge(this_STRING_3);
                    		

                    			newLeafNode(this_STRING_3, grammarAccess.getDotIDAccess().getSTRINGTerminalRuleCall_3());
                    		

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
    // $ANTLR end "ruleDotID"


    // $ANTLR start "ruleEdgeOperator"
    // InternalGraphvizDot.g:1037:1: ruleEdgeOperator returns [Enumerator current=null] : ( (enumLiteral_0= '->' ) | (enumLiteral_1= '--' ) ) ;
    public final Enumerator ruleEdgeOperator() throws RecognitionException {
        Enumerator current = null;

        Token enumLiteral_0=null;
        Token enumLiteral_1=null;


        	enterRule();

        try {
            // InternalGraphvizDot.g:1043:2: ( ( (enumLiteral_0= '->' ) | (enumLiteral_1= '--' ) ) )
            // InternalGraphvizDot.g:1044:2: ( (enumLiteral_0= '->' ) | (enumLiteral_1= '--' ) )
            {
            // InternalGraphvizDot.g:1044:2: ( (enumLiteral_0= '->' ) | (enumLiteral_1= '--' ) )
            int alt27=2;
            int LA27_0 = input.LA(1);

            if ( (LA27_0==23) ) {
                alt27=1;
            }
            else if ( (LA27_0==24) ) {
                alt27=2;
            }
            else {
                NoViableAltException nvae =
                    new NoViableAltException("", 27, 0, input);

                throw nvae;
            }
            switch (alt27) {
                case 1 :
                    // InternalGraphvizDot.g:1045:3: (enumLiteral_0= '->' )
                    {
                    // InternalGraphvizDot.g:1045:3: (enumLiteral_0= '->' )
                    // InternalGraphvizDot.g:1046:4: enumLiteral_0= '->'
                    {
                    enumLiteral_0=(Token)match(input,23,FOLLOW_2); 

                    				current = grammarAccess.getEdgeOperatorAccess().getDirectedEnumLiteralDeclaration_0().getEnumLiteral().getInstance();
                    				newLeafNode(enumLiteral_0, grammarAccess.getEdgeOperatorAccess().getDirectedEnumLiteralDeclaration_0());
                    			

                    }


                    }
                    break;
                case 2 :
                    // InternalGraphvizDot.g:1053:3: (enumLiteral_1= '--' )
                    {
                    // InternalGraphvizDot.g:1053:3: (enumLiteral_1= '--' )
                    // InternalGraphvizDot.g:1054:4: enumLiteral_1= '--'
                    {
                    enumLiteral_1=(Token)match(input,24,FOLLOW_2); 

                    				current = grammarAccess.getEdgeOperatorAccess().getUndirectedEnumLiteralDeclaration_1().getEnumLiteral().getInstance();
                    				newLeafNode(enumLiteral_1, grammarAccess.getEdgeOperatorAccess().getUndirectedEnumLiteralDeclaration_1());
                    			

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
    // $ANTLR end "ruleEdgeOperator"


    // $ANTLR start "ruleGraphType"
    // InternalGraphvizDot.g:1064:1: ruleGraphType returns [Enumerator current=null] : ( (enumLiteral_0= 'graph' ) | (enumLiteral_1= 'digraph' ) ) ;
    public final Enumerator ruleGraphType() throws RecognitionException {
        Enumerator current = null;

        Token enumLiteral_0=null;
        Token enumLiteral_1=null;


        	enterRule();

        try {
            // InternalGraphvizDot.g:1070:2: ( ( (enumLiteral_0= 'graph' ) | (enumLiteral_1= 'digraph' ) ) )
            // InternalGraphvizDot.g:1071:2: ( (enumLiteral_0= 'graph' ) | (enumLiteral_1= 'digraph' ) )
            {
            // InternalGraphvizDot.g:1071:2: ( (enumLiteral_0= 'graph' ) | (enumLiteral_1= 'digraph' ) )
            int alt28=2;
            int LA28_0 = input.LA(1);

            if ( (LA28_0==25) ) {
                alt28=1;
            }
            else if ( (LA28_0==26) ) {
                alt28=2;
            }
            else {
                NoViableAltException nvae =
                    new NoViableAltException("", 28, 0, input);

                throw nvae;
            }
            switch (alt28) {
                case 1 :
                    // InternalGraphvizDot.g:1072:3: (enumLiteral_0= 'graph' )
                    {
                    // InternalGraphvizDot.g:1072:3: (enumLiteral_0= 'graph' )
                    // InternalGraphvizDot.g:1073:4: enumLiteral_0= 'graph'
                    {
                    enumLiteral_0=(Token)match(input,25,FOLLOW_2); 

                    				current = grammarAccess.getGraphTypeAccess().getGraphEnumLiteralDeclaration_0().getEnumLiteral().getInstance();
                    				newLeafNode(enumLiteral_0, grammarAccess.getGraphTypeAccess().getGraphEnumLiteralDeclaration_0());
                    			

                    }


                    }
                    break;
                case 2 :
                    // InternalGraphvizDot.g:1080:3: (enumLiteral_1= 'digraph' )
                    {
                    // InternalGraphvizDot.g:1080:3: (enumLiteral_1= 'digraph' )
                    // InternalGraphvizDot.g:1081:4: enumLiteral_1= 'digraph'
                    {
                    enumLiteral_1=(Token)match(input,26,FOLLOW_2); 

                    				current = grammarAccess.getGraphTypeAccess().getDigraphEnumLiteralDeclaration_1().getEnumLiteral().getInstance();
                    				newLeafNode(enumLiteral_1, grammarAccess.getGraphTypeAccess().getDigraphEnumLiteralDeclaration_1());
                    			

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
    // $ANTLR end "ruleGraphType"


    // $ANTLR start "ruleAttributeType"
    // InternalGraphvizDot.g:1091:1: ruleAttributeType returns [Enumerator current=null] : ( (enumLiteral_0= 'graph' ) | (enumLiteral_1= 'node' ) | (enumLiteral_2= 'edge' ) ) ;
    public final Enumerator ruleAttributeType() throws RecognitionException {
        Enumerator current = null;

        Token enumLiteral_0=null;
        Token enumLiteral_1=null;
        Token enumLiteral_2=null;


        	enterRule();

        try {
            // InternalGraphvizDot.g:1097:2: ( ( (enumLiteral_0= 'graph' ) | (enumLiteral_1= 'node' ) | (enumLiteral_2= 'edge' ) ) )
            // InternalGraphvizDot.g:1098:2: ( (enumLiteral_0= 'graph' ) | (enumLiteral_1= 'node' ) | (enumLiteral_2= 'edge' ) )
            {
            // InternalGraphvizDot.g:1098:2: ( (enumLiteral_0= 'graph' ) | (enumLiteral_1= 'node' ) | (enumLiteral_2= 'edge' ) )
            int alt29=3;
            switch ( input.LA(1) ) {
            case 25:
                {
                alt29=1;
                }
                break;
            case 27:
                {
                alt29=2;
                }
                break;
            case 28:
                {
                alt29=3;
                }
                break;
            default:
                NoViableAltException nvae =
                    new NoViableAltException("", 29, 0, input);

                throw nvae;
            }

            switch (alt29) {
                case 1 :
                    // InternalGraphvizDot.g:1099:3: (enumLiteral_0= 'graph' )
                    {
                    // InternalGraphvizDot.g:1099:3: (enumLiteral_0= 'graph' )
                    // InternalGraphvizDot.g:1100:4: enumLiteral_0= 'graph'
                    {
                    enumLiteral_0=(Token)match(input,25,FOLLOW_2); 

                    				current = grammarAccess.getAttributeTypeAccess().getGraphEnumLiteralDeclaration_0().getEnumLiteral().getInstance();
                    				newLeafNode(enumLiteral_0, grammarAccess.getAttributeTypeAccess().getGraphEnumLiteralDeclaration_0());
                    			

                    }


                    }
                    break;
                case 2 :
                    // InternalGraphvizDot.g:1107:3: (enumLiteral_1= 'node' )
                    {
                    // InternalGraphvizDot.g:1107:3: (enumLiteral_1= 'node' )
                    // InternalGraphvizDot.g:1108:4: enumLiteral_1= 'node'
                    {
                    enumLiteral_1=(Token)match(input,27,FOLLOW_2); 

                    				current = grammarAccess.getAttributeTypeAccess().getNodeEnumLiteralDeclaration_1().getEnumLiteral().getInstance();
                    				newLeafNode(enumLiteral_1, grammarAccess.getAttributeTypeAccess().getNodeEnumLiteralDeclaration_1());
                    			

                    }


                    }
                    break;
                case 3 :
                    // InternalGraphvizDot.g:1115:3: (enumLiteral_2= 'edge' )
                    {
                    // InternalGraphvizDot.g:1115:3: (enumLiteral_2= 'edge' )
                    // InternalGraphvizDot.g:1116:4: enumLiteral_2= 'edge'
                    {
                    enumLiteral_2=(Token)match(input,28,FOLLOW_2); 

                    				current = grammarAccess.getAttributeTypeAccess().getEdgeEnumLiteralDeclaration_2().getEnumLiteral().getInstance();
                    				newLeafNode(enumLiteral_2, grammarAccess.getAttributeTypeAccess().getEdgeEnumLiteralDeclaration_2());
                    			

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
    // $ANTLR end "ruleAttributeType"

    // Delegated rules


    protected DFA5 dfa5 = new DFA5(this);
    static final String dfa_1s = "\21\uffff";
    static final String dfa_2s = "\1\uffff\4\12\6\uffff\4\12\1\uffff\1\12";
    static final String dfa_3s = "\5\4\3\uffff\1\4\2\uffff\6\4";
    static final String dfa_4s = "\5\34\3\uffff\1\7\2\uffff\4\34\1\4\1\34";
    static final String dfa_5s = "\5\uffff\1\3\1\5\1\4\1\uffff\1\2\1\1\6\uffff";
    static final String dfa_6s = "\21\uffff}>";
    static final String[] dfa_7s = {
            "\1\1\1\2\1\3\1\4\6\uffff\1\6\7\uffff\1\6\2\uffff\1\5\1\uffff\2\5",
            "\4\12\6\uffff\3\12\1\7\1\12\2\uffff\1\10\1\12\2\11\1\12\1\uffff\2\12",
            "\4\12\6\uffff\3\12\1\7\1\12\2\uffff\1\10\1\12\2\11\1\12\1\uffff\2\12",
            "\4\12\6\uffff\3\12\1\7\1\12\2\uffff\1\10\1\12\2\11\1\12\1\uffff\2\12",
            "\4\12\6\uffff\3\12\1\7\1\12\2\uffff\1\10\1\12\2\11\1\12\1\uffff\2\12",
            "",
            "",
            "",
            "\1\13\1\14\1\15\1\16",
            "",
            "",
            "\4\12\6\uffff\3\12\1\uffff\1\12\2\uffff\1\17\1\12\2\11\1\12\1\uffff\2\12",
            "\4\12\6\uffff\3\12\1\uffff\1\12\2\uffff\1\17\1\12\2\11\1\12\1\uffff\2\12",
            "\4\12\6\uffff\3\12\1\uffff\1\12\2\uffff\1\17\1\12\2\11\1\12\1\uffff\2\12",
            "\4\12\6\uffff\3\12\1\uffff\1\12\2\uffff\1\17\1\12\2\11\1\12\1\uffff\2\12",
            "\1\20",
            "\4\12\6\uffff\3\12\1\uffff\1\12\3\uffff\1\12\2\11\1\12\1\uffff\2\12"
    };

    static final short[] dfa_1 = DFA.unpackEncodedString(dfa_1s);
    static final short[] dfa_2 = DFA.unpackEncodedString(dfa_2s);
    static final char[] dfa_3 = DFA.unpackEncodedStringToUnsignedChars(dfa_3s);
    static final char[] dfa_4 = DFA.unpackEncodedStringToUnsignedChars(dfa_4s);
    static final short[] dfa_5 = DFA.unpackEncodedString(dfa_5s);
    static final short[] dfa_6 = DFA.unpackEncodedString(dfa_6s);
    static final short[][] dfa_7 = unpackEncodedStringArray(dfa_7s);

    class DFA5 extends DFA {

        public DFA5(BaseRecognizer recognizer) {
            this.recognizer = recognizer;
            this.decisionNumber = 5;
            this.eot = dfa_1;
            this.eof = dfa_2;
            this.min = dfa_3;
            this.max = dfa_4;
            this.accept = dfa_5;
            this.special = dfa_6;
            this.transition = dfa_7;
        }
        public String getDescription() {
            return "221:3: (this_NodeStatement_0= ruleNodeStatement | this_EdgeStatement_1= ruleEdgeStatement | this_AttributeStatement_2= ruleAttributeStatement | this_Attribute_3= ruleAttribute | this_Subgraph_4= ruleSubgraph )";
        }
    }
 

    public static final BitSet FOLLOW_1 = new BitSet(new long[]{0x0000000000000000L});
    public static final BitSet FOLLOW_2 = new BitSet(new long[]{0x0000000000000002L});
    public static final BitSet FOLLOW_3 = new BitSet(new long[]{0x0000000006002002L});
    public static final BitSet FOLLOW_4 = new BitSet(new long[]{0x00000000060060F0L});
    public static final BitSet FOLLOW_5 = new BitSet(new long[]{0x00000000000040F0L});
    public static final BitSet FOLLOW_6 = new BitSet(new long[]{0x0000000000004000L});
    public static final BitSet FOLLOW_7 = new BitSet(new long[]{0x000000001A40C0F0L});
    public static final BitSet FOLLOW_8 = new BitSet(new long[]{0x0000000000010002L});
    public static final BitSet FOLLOW_9 = new BitSet(new long[]{0x0000000000020000L});
    public static final BitSet FOLLOW_10 = new BitSet(new long[]{0x00000000000000F0L});
    public static final BitSet FOLLOW_11 = new BitSet(new long[]{0x0000000000040002L});
    public static final BitSet FOLLOW_12 = new BitSet(new long[]{0x00000000001000F0L});
    public static final BitSet FOLLOW_13 = new BitSet(new long[]{0x00000000001800F0L});
    public static final BitSet FOLLOW_14 = new BitSet(new long[]{0x0000000000200002L});
    public static final BitSet FOLLOW_15 = new BitSet(new long[]{0x0000000001800000L});
    public static final BitSet FOLLOW_16 = new BitSet(new long[]{0x0000000001840002L});
    public static final BitSet FOLLOW_17 = new BitSet(new long[]{0x000000001A4040F0L});
    public static final BitSet FOLLOW_18 = new BitSet(new long[]{0x0000000000040000L});
    public static final BitSet FOLLOW_19 = new BitSet(new long[]{0x0000000000004010L});
    public static final BitSet FOLLOW_20 = new BitSet(new long[]{0x0000000000020002L});
    public static final BitSet FOLLOW_21 = new BitSet(new long[]{0x0000000000000010L});

}