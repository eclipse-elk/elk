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
package org.eclipse.elk.graph.text.ide.contentassist.antlr.internal;

import java.io.InputStream;
import org.eclipse.xtext.*;
import org.eclipse.xtext.parser.*;
import org.eclipse.xtext.parser.impl.*;
import org.eclipse.emf.ecore.util.EcoreUtil;
import org.eclipse.emf.ecore.EObject;
import org.eclipse.xtext.parser.antlr.XtextTokenStream;
import org.eclipse.xtext.parser.antlr.XtextTokenStream.HiddenTokens;
import org.eclipse.xtext.ide.editor.contentassist.antlr.internal.AbstractInternalContentAssistParser;
import org.eclipse.xtext.ide.editor.contentassist.antlr.internal.DFA;
import org.eclipse.elk.graph.text.services.ElkGraphGrammarAccess;



import org.antlr.runtime.*;
import java.util.Stack;
import java.util.List;
import java.util.ArrayList;

/*******************************************************************************
 * Copyright (c) 2016 TypeFox GmbH (http://www.typefox.io) and others.
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
public class InternalElkGraphParser extends AbstractInternalContentAssistParser {
    public static final String[] tokenNames = new String[] {
        "<invalid>", "<EOR>", "<DOWN>", "<UP>", "RULE_STRING", "RULE_SIGNED_INT", "RULE_FLOAT", "RULE_ID", "RULE_INT", "RULE_ML_COMMENT", "RULE_SL_COMMENT", "RULE_WS", "RULE_ANY_OTHER", "'null'", "'true'", "'false'", "'graph'", "'node'", "'{'", "'}'", "'label'", "':'", "'port'", "'layout'", "'['", "']'", "'position'", "','", "'size'", "'edge'", "'->'", "'incoming'", "'outgoing'", "'start'", "'end'", "'bends'", "'|'", "'section'", "'.'"
    };
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


        public InternalElkGraphParser(TokenStream input) {
            this(input, new RecognizerSharedState());
        }
        public InternalElkGraphParser(TokenStream input, RecognizerSharedState state) {
            super(input, state);
             
        }
        

    public String[] getTokenNames() { return InternalElkGraphParser.tokenNames; }
    public String getGrammarFileName() { return "InternalElkGraph.g"; }


    	private ElkGraphGrammarAccess grammarAccess;

    	public void setGrammarAccess(ElkGraphGrammarAccess grammarAccess) {
    		this.grammarAccess = grammarAccess;
    	}

    	@Override
    	protected Grammar getGrammar() {
    		return grammarAccess.getGrammar();
    	}

    	@Override
    	protected String getValueForTokenName(String tokenName) {
    		return tokenName;
    	}



    // $ANTLR start "entryRuleRootNode"
    // InternalElkGraph.g:57:1: entryRuleRootNode : ruleRootNode EOF ;
    public final void entryRuleRootNode() throws RecognitionException {
        try {
            // InternalElkGraph.g:58:1: ( ruleRootNode EOF )
            // InternalElkGraph.g:59:1: ruleRootNode EOF
            {
             before(grammarAccess.getRootNodeRule()); 
            pushFollow(FOLLOW_1);
            ruleRootNode();

            state._fsp--;

             after(grammarAccess.getRootNodeRule()); 
            match(input,EOF,FOLLOW_2); 

            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {
        }
        return ;
    }
    // $ANTLR end "entryRuleRootNode"


    // $ANTLR start "ruleRootNode"
    // InternalElkGraph.g:66:1: ruleRootNode : ( ( rule__RootNode__Group__0 ) ) ;
    public final void ruleRootNode() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:70:2: ( ( ( rule__RootNode__Group__0 ) ) )
            // InternalElkGraph.g:71:2: ( ( rule__RootNode__Group__0 ) )
            {
            // InternalElkGraph.g:71:2: ( ( rule__RootNode__Group__0 ) )
            // InternalElkGraph.g:72:3: ( rule__RootNode__Group__0 )
            {
             before(grammarAccess.getRootNodeAccess().getGroup()); 
            // InternalElkGraph.g:73:3: ( rule__RootNode__Group__0 )
            // InternalElkGraph.g:73:4: rule__RootNode__Group__0
            {
            pushFollow(FOLLOW_2);
            rule__RootNode__Group__0();

            state._fsp--;


            }

             after(grammarAccess.getRootNodeAccess().getGroup()); 

            }


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "ruleRootNode"


    // $ANTLR start "entryRuleElkNode"
    // InternalElkGraph.g:82:1: entryRuleElkNode : ruleElkNode EOF ;
    public final void entryRuleElkNode() throws RecognitionException {
        try {
            // InternalElkGraph.g:83:1: ( ruleElkNode EOF )
            // InternalElkGraph.g:84:1: ruleElkNode EOF
            {
             before(grammarAccess.getElkNodeRule()); 
            pushFollow(FOLLOW_1);
            ruleElkNode();

            state._fsp--;

             after(grammarAccess.getElkNodeRule()); 
            match(input,EOF,FOLLOW_2); 

            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {
        }
        return ;
    }
    // $ANTLR end "entryRuleElkNode"


    // $ANTLR start "ruleElkNode"
    // InternalElkGraph.g:91:1: ruleElkNode : ( ( rule__ElkNode__Group__0 ) ) ;
    public final void ruleElkNode() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:95:2: ( ( ( rule__ElkNode__Group__0 ) ) )
            // InternalElkGraph.g:96:2: ( ( rule__ElkNode__Group__0 ) )
            {
            // InternalElkGraph.g:96:2: ( ( rule__ElkNode__Group__0 ) )
            // InternalElkGraph.g:97:3: ( rule__ElkNode__Group__0 )
            {
             before(grammarAccess.getElkNodeAccess().getGroup()); 
            // InternalElkGraph.g:98:3: ( rule__ElkNode__Group__0 )
            // InternalElkGraph.g:98:4: rule__ElkNode__Group__0
            {
            pushFollow(FOLLOW_2);
            rule__ElkNode__Group__0();

            state._fsp--;


            }

             after(grammarAccess.getElkNodeAccess().getGroup()); 

            }


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "ruleElkNode"


    // $ANTLR start "entryRuleElkLabel"
    // InternalElkGraph.g:107:1: entryRuleElkLabel : ruleElkLabel EOF ;
    public final void entryRuleElkLabel() throws RecognitionException {
        try {
            // InternalElkGraph.g:108:1: ( ruleElkLabel EOF )
            // InternalElkGraph.g:109:1: ruleElkLabel EOF
            {
             before(grammarAccess.getElkLabelRule()); 
            pushFollow(FOLLOW_1);
            ruleElkLabel();

            state._fsp--;

             after(grammarAccess.getElkLabelRule()); 
            match(input,EOF,FOLLOW_2); 

            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {
        }
        return ;
    }
    // $ANTLR end "entryRuleElkLabel"


    // $ANTLR start "ruleElkLabel"
    // InternalElkGraph.g:116:1: ruleElkLabel : ( ( rule__ElkLabel__Group__0 ) ) ;
    public final void ruleElkLabel() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:120:2: ( ( ( rule__ElkLabel__Group__0 ) ) )
            // InternalElkGraph.g:121:2: ( ( rule__ElkLabel__Group__0 ) )
            {
            // InternalElkGraph.g:121:2: ( ( rule__ElkLabel__Group__0 ) )
            // InternalElkGraph.g:122:3: ( rule__ElkLabel__Group__0 )
            {
             before(grammarAccess.getElkLabelAccess().getGroup()); 
            // InternalElkGraph.g:123:3: ( rule__ElkLabel__Group__0 )
            // InternalElkGraph.g:123:4: rule__ElkLabel__Group__0
            {
            pushFollow(FOLLOW_2);
            rule__ElkLabel__Group__0();

            state._fsp--;


            }

             after(grammarAccess.getElkLabelAccess().getGroup()); 

            }


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "ruleElkLabel"


    // $ANTLR start "entryRuleElkPort"
    // InternalElkGraph.g:132:1: entryRuleElkPort : ruleElkPort EOF ;
    public final void entryRuleElkPort() throws RecognitionException {
        try {
            // InternalElkGraph.g:133:1: ( ruleElkPort EOF )
            // InternalElkGraph.g:134:1: ruleElkPort EOF
            {
             before(grammarAccess.getElkPortRule()); 
            pushFollow(FOLLOW_1);
            ruleElkPort();

            state._fsp--;

             after(grammarAccess.getElkPortRule()); 
            match(input,EOF,FOLLOW_2); 

            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {
        }
        return ;
    }
    // $ANTLR end "entryRuleElkPort"


    // $ANTLR start "ruleElkPort"
    // InternalElkGraph.g:141:1: ruleElkPort : ( ( rule__ElkPort__Group__0 ) ) ;
    public final void ruleElkPort() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:145:2: ( ( ( rule__ElkPort__Group__0 ) ) )
            // InternalElkGraph.g:146:2: ( ( rule__ElkPort__Group__0 ) )
            {
            // InternalElkGraph.g:146:2: ( ( rule__ElkPort__Group__0 ) )
            // InternalElkGraph.g:147:3: ( rule__ElkPort__Group__0 )
            {
             before(grammarAccess.getElkPortAccess().getGroup()); 
            // InternalElkGraph.g:148:3: ( rule__ElkPort__Group__0 )
            // InternalElkGraph.g:148:4: rule__ElkPort__Group__0
            {
            pushFollow(FOLLOW_2);
            rule__ElkPort__Group__0();

            state._fsp--;


            }

             after(grammarAccess.getElkPortAccess().getGroup()); 

            }


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "ruleElkPort"


    // $ANTLR start "ruleShapeLayout"
    // InternalElkGraph.g:158:1: ruleShapeLayout : ( ( rule__ShapeLayout__Group__0 ) ) ;
    public final void ruleShapeLayout() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:162:2: ( ( ( rule__ShapeLayout__Group__0 ) ) )
            // InternalElkGraph.g:163:2: ( ( rule__ShapeLayout__Group__0 ) )
            {
            // InternalElkGraph.g:163:2: ( ( rule__ShapeLayout__Group__0 ) )
            // InternalElkGraph.g:164:3: ( rule__ShapeLayout__Group__0 )
            {
             before(grammarAccess.getShapeLayoutAccess().getGroup()); 
            // InternalElkGraph.g:165:3: ( rule__ShapeLayout__Group__0 )
            // InternalElkGraph.g:165:4: rule__ShapeLayout__Group__0
            {
            pushFollow(FOLLOW_2);
            rule__ShapeLayout__Group__0();

            state._fsp--;


            }

             after(grammarAccess.getShapeLayoutAccess().getGroup()); 

            }


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "ruleShapeLayout"


    // $ANTLR start "entryRuleElkEdge"
    // InternalElkGraph.g:174:1: entryRuleElkEdge : ruleElkEdge EOF ;
    public final void entryRuleElkEdge() throws RecognitionException {
        try {
            // InternalElkGraph.g:175:1: ( ruleElkEdge EOF )
            // InternalElkGraph.g:176:1: ruleElkEdge EOF
            {
             before(grammarAccess.getElkEdgeRule()); 
            pushFollow(FOLLOW_1);
            ruleElkEdge();

            state._fsp--;

             after(grammarAccess.getElkEdgeRule()); 
            match(input,EOF,FOLLOW_2); 

            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {
        }
        return ;
    }
    // $ANTLR end "entryRuleElkEdge"


    // $ANTLR start "ruleElkEdge"
    // InternalElkGraph.g:183:1: ruleElkEdge : ( ( rule__ElkEdge__Group__0 ) ) ;
    public final void ruleElkEdge() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:187:2: ( ( ( rule__ElkEdge__Group__0 ) ) )
            // InternalElkGraph.g:188:2: ( ( rule__ElkEdge__Group__0 ) )
            {
            // InternalElkGraph.g:188:2: ( ( rule__ElkEdge__Group__0 ) )
            // InternalElkGraph.g:189:3: ( rule__ElkEdge__Group__0 )
            {
             before(grammarAccess.getElkEdgeAccess().getGroup()); 
            // InternalElkGraph.g:190:3: ( rule__ElkEdge__Group__0 )
            // InternalElkGraph.g:190:4: rule__ElkEdge__Group__0
            {
            pushFollow(FOLLOW_2);
            rule__ElkEdge__Group__0();

            state._fsp--;


            }

             after(grammarAccess.getElkEdgeAccess().getGroup()); 

            }


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "ruleElkEdge"


    // $ANTLR start "ruleEdgeLayout"
    // InternalElkGraph.g:200:1: ruleEdgeLayout : ( ( rule__EdgeLayout__Group__0 ) ) ;
    public final void ruleEdgeLayout() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:204:2: ( ( ( rule__EdgeLayout__Group__0 ) ) )
            // InternalElkGraph.g:205:2: ( ( rule__EdgeLayout__Group__0 ) )
            {
            // InternalElkGraph.g:205:2: ( ( rule__EdgeLayout__Group__0 ) )
            // InternalElkGraph.g:206:3: ( rule__EdgeLayout__Group__0 )
            {
             before(grammarAccess.getEdgeLayoutAccess().getGroup()); 
            // InternalElkGraph.g:207:3: ( rule__EdgeLayout__Group__0 )
            // InternalElkGraph.g:207:4: rule__EdgeLayout__Group__0
            {
            pushFollow(FOLLOW_2);
            rule__EdgeLayout__Group__0();

            state._fsp--;


            }

             after(grammarAccess.getEdgeLayoutAccess().getGroup()); 

            }


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "ruleEdgeLayout"


    // $ANTLR start "entryRuleElkSingleEdgeSection"
    // InternalElkGraph.g:216:1: entryRuleElkSingleEdgeSection : ruleElkSingleEdgeSection EOF ;
    public final void entryRuleElkSingleEdgeSection() throws RecognitionException {
        try {
            // InternalElkGraph.g:217:1: ( ruleElkSingleEdgeSection EOF )
            // InternalElkGraph.g:218:1: ruleElkSingleEdgeSection EOF
            {
             before(grammarAccess.getElkSingleEdgeSectionRule()); 
            pushFollow(FOLLOW_1);
            ruleElkSingleEdgeSection();

            state._fsp--;

             after(grammarAccess.getElkSingleEdgeSectionRule()); 
            match(input,EOF,FOLLOW_2); 

            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {
        }
        return ;
    }
    // $ANTLR end "entryRuleElkSingleEdgeSection"


    // $ANTLR start "ruleElkSingleEdgeSection"
    // InternalElkGraph.g:225:1: ruleElkSingleEdgeSection : ( ( rule__ElkSingleEdgeSection__Group__0 ) ) ;
    public final void ruleElkSingleEdgeSection() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:229:2: ( ( ( rule__ElkSingleEdgeSection__Group__0 ) ) )
            // InternalElkGraph.g:230:2: ( ( rule__ElkSingleEdgeSection__Group__0 ) )
            {
            // InternalElkGraph.g:230:2: ( ( rule__ElkSingleEdgeSection__Group__0 ) )
            // InternalElkGraph.g:231:3: ( rule__ElkSingleEdgeSection__Group__0 )
            {
             before(grammarAccess.getElkSingleEdgeSectionAccess().getGroup()); 
            // InternalElkGraph.g:232:3: ( rule__ElkSingleEdgeSection__Group__0 )
            // InternalElkGraph.g:232:4: rule__ElkSingleEdgeSection__Group__0
            {
            pushFollow(FOLLOW_2);
            rule__ElkSingleEdgeSection__Group__0();

            state._fsp--;


            }

             after(grammarAccess.getElkSingleEdgeSectionAccess().getGroup()); 

            }


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "ruleElkSingleEdgeSection"


    // $ANTLR start "entryRuleElkEdgeSection"
    // InternalElkGraph.g:241:1: entryRuleElkEdgeSection : ruleElkEdgeSection EOF ;
    public final void entryRuleElkEdgeSection() throws RecognitionException {
        try {
            // InternalElkGraph.g:242:1: ( ruleElkEdgeSection EOF )
            // InternalElkGraph.g:243:1: ruleElkEdgeSection EOF
            {
             before(grammarAccess.getElkEdgeSectionRule()); 
            pushFollow(FOLLOW_1);
            ruleElkEdgeSection();

            state._fsp--;

             after(grammarAccess.getElkEdgeSectionRule()); 
            match(input,EOF,FOLLOW_2); 

            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {
        }
        return ;
    }
    // $ANTLR end "entryRuleElkEdgeSection"


    // $ANTLR start "ruleElkEdgeSection"
    // InternalElkGraph.g:250:1: ruleElkEdgeSection : ( ( rule__ElkEdgeSection__Group__0 ) ) ;
    public final void ruleElkEdgeSection() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:254:2: ( ( ( rule__ElkEdgeSection__Group__0 ) ) )
            // InternalElkGraph.g:255:2: ( ( rule__ElkEdgeSection__Group__0 ) )
            {
            // InternalElkGraph.g:255:2: ( ( rule__ElkEdgeSection__Group__0 ) )
            // InternalElkGraph.g:256:3: ( rule__ElkEdgeSection__Group__0 )
            {
             before(grammarAccess.getElkEdgeSectionAccess().getGroup()); 
            // InternalElkGraph.g:257:3: ( rule__ElkEdgeSection__Group__0 )
            // InternalElkGraph.g:257:4: rule__ElkEdgeSection__Group__0
            {
            pushFollow(FOLLOW_2);
            rule__ElkEdgeSection__Group__0();

            state._fsp--;


            }

             after(grammarAccess.getElkEdgeSectionAccess().getGroup()); 

            }


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "ruleElkEdgeSection"


    // $ANTLR start "entryRuleElkBendPoint"
    // InternalElkGraph.g:266:1: entryRuleElkBendPoint : ruleElkBendPoint EOF ;
    public final void entryRuleElkBendPoint() throws RecognitionException {
        try {
            // InternalElkGraph.g:267:1: ( ruleElkBendPoint EOF )
            // InternalElkGraph.g:268:1: ruleElkBendPoint EOF
            {
             before(grammarAccess.getElkBendPointRule()); 
            pushFollow(FOLLOW_1);
            ruleElkBendPoint();

            state._fsp--;

             after(grammarAccess.getElkBendPointRule()); 
            match(input,EOF,FOLLOW_2); 

            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {
        }
        return ;
    }
    // $ANTLR end "entryRuleElkBendPoint"


    // $ANTLR start "ruleElkBendPoint"
    // InternalElkGraph.g:275:1: ruleElkBendPoint : ( ( rule__ElkBendPoint__Group__0 ) ) ;
    public final void ruleElkBendPoint() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:279:2: ( ( ( rule__ElkBendPoint__Group__0 ) ) )
            // InternalElkGraph.g:280:2: ( ( rule__ElkBendPoint__Group__0 ) )
            {
            // InternalElkGraph.g:280:2: ( ( rule__ElkBendPoint__Group__0 ) )
            // InternalElkGraph.g:281:3: ( rule__ElkBendPoint__Group__0 )
            {
             before(grammarAccess.getElkBendPointAccess().getGroup()); 
            // InternalElkGraph.g:282:3: ( rule__ElkBendPoint__Group__0 )
            // InternalElkGraph.g:282:4: rule__ElkBendPoint__Group__0
            {
            pushFollow(FOLLOW_2);
            rule__ElkBendPoint__Group__0();

            state._fsp--;


            }

             after(grammarAccess.getElkBendPointAccess().getGroup()); 

            }


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "ruleElkBendPoint"


    // $ANTLR start "entryRuleQualifiedId"
    // InternalElkGraph.g:291:1: entryRuleQualifiedId : ruleQualifiedId EOF ;
    public final void entryRuleQualifiedId() throws RecognitionException {
        try {
            // InternalElkGraph.g:292:1: ( ruleQualifiedId EOF )
            // InternalElkGraph.g:293:1: ruleQualifiedId EOF
            {
             before(grammarAccess.getQualifiedIdRule()); 
            pushFollow(FOLLOW_1);
            ruleQualifiedId();

            state._fsp--;

             after(grammarAccess.getQualifiedIdRule()); 
            match(input,EOF,FOLLOW_2); 

            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {
        }
        return ;
    }
    // $ANTLR end "entryRuleQualifiedId"


    // $ANTLR start "ruleQualifiedId"
    // InternalElkGraph.g:300:1: ruleQualifiedId : ( ( rule__QualifiedId__Group__0 ) ) ;
    public final void ruleQualifiedId() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:304:2: ( ( ( rule__QualifiedId__Group__0 ) ) )
            // InternalElkGraph.g:305:2: ( ( rule__QualifiedId__Group__0 ) )
            {
            // InternalElkGraph.g:305:2: ( ( rule__QualifiedId__Group__0 ) )
            // InternalElkGraph.g:306:3: ( rule__QualifiedId__Group__0 )
            {
             before(grammarAccess.getQualifiedIdAccess().getGroup()); 
            // InternalElkGraph.g:307:3: ( rule__QualifiedId__Group__0 )
            // InternalElkGraph.g:307:4: rule__QualifiedId__Group__0
            {
            pushFollow(FOLLOW_2);
            rule__QualifiedId__Group__0();

            state._fsp--;


            }

             after(grammarAccess.getQualifiedIdAccess().getGroup()); 

            }


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "ruleQualifiedId"


    // $ANTLR start "entryRuleNumber"
    // InternalElkGraph.g:316:1: entryRuleNumber : ruleNumber EOF ;
    public final void entryRuleNumber() throws RecognitionException {
        try {
            // InternalElkGraph.g:317:1: ( ruleNumber EOF )
            // InternalElkGraph.g:318:1: ruleNumber EOF
            {
             before(grammarAccess.getNumberRule()); 
            pushFollow(FOLLOW_1);
            ruleNumber();

            state._fsp--;

             after(grammarAccess.getNumberRule()); 
            match(input,EOF,FOLLOW_2); 

            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {
        }
        return ;
    }
    // $ANTLR end "entryRuleNumber"


    // $ANTLR start "ruleNumber"
    // InternalElkGraph.g:325:1: ruleNumber : ( ( rule__Number__Alternatives ) ) ;
    public final void ruleNumber() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:329:2: ( ( ( rule__Number__Alternatives ) ) )
            // InternalElkGraph.g:330:2: ( ( rule__Number__Alternatives ) )
            {
            // InternalElkGraph.g:330:2: ( ( rule__Number__Alternatives ) )
            // InternalElkGraph.g:331:3: ( rule__Number__Alternatives )
            {
             before(grammarAccess.getNumberAccess().getAlternatives()); 
            // InternalElkGraph.g:332:3: ( rule__Number__Alternatives )
            // InternalElkGraph.g:332:4: rule__Number__Alternatives
            {
            pushFollow(FOLLOW_2);
            rule__Number__Alternatives();

            state._fsp--;


            }

             after(grammarAccess.getNumberAccess().getAlternatives()); 

            }


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "ruleNumber"


    // $ANTLR start "entryRuleProperty"
    // InternalElkGraph.g:341:1: entryRuleProperty : ruleProperty EOF ;
    public final void entryRuleProperty() throws RecognitionException {
        try {
            // InternalElkGraph.g:342:1: ( ruleProperty EOF )
            // InternalElkGraph.g:343:1: ruleProperty EOF
            {
             before(grammarAccess.getPropertyRule()); 
            pushFollow(FOLLOW_1);
            ruleProperty();

            state._fsp--;

             after(grammarAccess.getPropertyRule()); 
            match(input,EOF,FOLLOW_2); 

            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {
        }
        return ;
    }
    // $ANTLR end "entryRuleProperty"


    // $ANTLR start "ruleProperty"
    // InternalElkGraph.g:350:1: ruleProperty : ( ( rule__Property__Group__0 ) ) ;
    public final void ruleProperty() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:354:2: ( ( ( rule__Property__Group__0 ) ) )
            // InternalElkGraph.g:355:2: ( ( rule__Property__Group__0 ) )
            {
            // InternalElkGraph.g:355:2: ( ( rule__Property__Group__0 ) )
            // InternalElkGraph.g:356:3: ( rule__Property__Group__0 )
            {
             before(grammarAccess.getPropertyAccess().getGroup()); 
            // InternalElkGraph.g:357:3: ( rule__Property__Group__0 )
            // InternalElkGraph.g:357:4: rule__Property__Group__0
            {
            pushFollow(FOLLOW_2);
            rule__Property__Group__0();

            state._fsp--;


            }

             after(grammarAccess.getPropertyAccess().getGroup()); 

            }


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "ruleProperty"


    // $ANTLR start "entryRulePropertyKey"
    // InternalElkGraph.g:366:1: entryRulePropertyKey : rulePropertyKey EOF ;
    public final void entryRulePropertyKey() throws RecognitionException {
         
        	HiddenTokens myHiddenTokenState = ((XtextTokenStream)input).setHiddenTokens();

        try {
            // InternalElkGraph.g:370:1: ( rulePropertyKey EOF )
            // InternalElkGraph.g:371:1: rulePropertyKey EOF
            {
             before(grammarAccess.getPropertyKeyRule()); 
            pushFollow(FOLLOW_1);
            rulePropertyKey();

            state._fsp--;

             after(grammarAccess.getPropertyKeyRule()); 
            match(input,EOF,FOLLOW_2); 

            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	myHiddenTokenState.restore();

        }
        return ;
    }
    // $ANTLR end "entryRulePropertyKey"


    // $ANTLR start "rulePropertyKey"
    // InternalElkGraph.g:381:1: rulePropertyKey : ( ( rule__PropertyKey__Group__0 ) ) ;
    public final void rulePropertyKey() throws RecognitionException {

        		HiddenTokens myHiddenTokenState = ((XtextTokenStream)input).setHiddenTokens();
        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:386:2: ( ( ( rule__PropertyKey__Group__0 ) ) )
            // InternalElkGraph.g:387:2: ( ( rule__PropertyKey__Group__0 ) )
            {
            // InternalElkGraph.g:387:2: ( ( rule__PropertyKey__Group__0 ) )
            // InternalElkGraph.g:388:3: ( rule__PropertyKey__Group__0 )
            {
             before(grammarAccess.getPropertyKeyAccess().getGroup()); 
            // InternalElkGraph.g:389:3: ( rule__PropertyKey__Group__0 )
            // InternalElkGraph.g:389:4: rule__PropertyKey__Group__0
            {
            pushFollow(FOLLOW_2);
            rule__PropertyKey__Group__0();

            state._fsp--;


            }

             after(grammarAccess.getPropertyKeyAccess().getGroup()); 

            }


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);
            	myHiddenTokenState.restore();

        }
        return ;
    }
    // $ANTLR end "rulePropertyKey"


    // $ANTLR start "entryRuleStringValue"
    // InternalElkGraph.g:399:1: entryRuleStringValue : ruleStringValue EOF ;
    public final void entryRuleStringValue() throws RecognitionException {
        try {
            // InternalElkGraph.g:400:1: ( ruleStringValue EOF )
            // InternalElkGraph.g:401:1: ruleStringValue EOF
            {
             before(grammarAccess.getStringValueRule()); 
            pushFollow(FOLLOW_1);
            ruleStringValue();

            state._fsp--;

             after(grammarAccess.getStringValueRule()); 
            match(input,EOF,FOLLOW_2); 

            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {
        }
        return ;
    }
    // $ANTLR end "entryRuleStringValue"


    // $ANTLR start "ruleStringValue"
    // InternalElkGraph.g:408:1: ruleStringValue : ( RULE_STRING ) ;
    public final void ruleStringValue() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:412:2: ( ( RULE_STRING ) )
            // InternalElkGraph.g:413:2: ( RULE_STRING )
            {
            // InternalElkGraph.g:413:2: ( RULE_STRING )
            // InternalElkGraph.g:414:3: RULE_STRING
            {
             before(grammarAccess.getStringValueAccess().getSTRINGTerminalRuleCall()); 
            match(input,RULE_STRING,FOLLOW_2); 
             after(grammarAccess.getStringValueAccess().getSTRINGTerminalRuleCall()); 

            }


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "ruleStringValue"


    // $ANTLR start "entryRuleQualifiedIdValue"
    // InternalElkGraph.g:424:1: entryRuleQualifiedIdValue : ruleQualifiedIdValue EOF ;
    public final void entryRuleQualifiedIdValue() throws RecognitionException {
        try {
            // InternalElkGraph.g:425:1: ( ruleQualifiedIdValue EOF )
            // InternalElkGraph.g:426:1: ruleQualifiedIdValue EOF
            {
             before(grammarAccess.getQualifiedIdValueRule()); 
            pushFollow(FOLLOW_1);
            ruleQualifiedIdValue();

            state._fsp--;

             after(grammarAccess.getQualifiedIdValueRule()); 
            match(input,EOF,FOLLOW_2); 

            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {
        }
        return ;
    }
    // $ANTLR end "entryRuleQualifiedIdValue"


    // $ANTLR start "ruleQualifiedIdValue"
    // InternalElkGraph.g:433:1: ruleQualifiedIdValue : ( ruleQualifiedId ) ;
    public final void ruleQualifiedIdValue() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:437:2: ( ( ruleQualifiedId ) )
            // InternalElkGraph.g:438:2: ( ruleQualifiedId )
            {
            // InternalElkGraph.g:438:2: ( ruleQualifiedId )
            // InternalElkGraph.g:439:3: ruleQualifiedId
            {
             before(grammarAccess.getQualifiedIdValueAccess().getQualifiedIdParserRuleCall()); 
            pushFollow(FOLLOW_2);
            ruleQualifiedId();

            state._fsp--;

             after(grammarAccess.getQualifiedIdValueAccess().getQualifiedIdParserRuleCall()); 

            }


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "ruleQualifiedIdValue"


    // $ANTLR start "entryRuleNumberValue"
    // InternalElkGraph.g:449:1: entryRuleNumberValue : ruleNumberValue EOF ;
    public final void entryRuleNumberValue() throws RecognitionException {
        try {
            // InternalElkGraph.g:450:1: ( ruleNumberValue EOF )
            // InternalElkGraph.g:451:1: ruleNumberValue EOF
            {
             before(grammarAccess.getNumberValueRule()); 
            pushFollow(FOLLOW_1);
            ruleNumberValue();

            state._fsp--;

             after(grammarAccess.getNumberValueRule()); 
            match(input,EOF,FOLLOW_2); 

            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {
        }
        return ;
    }
    // $ANTLR end "entryRuleNumberValue"


    // $ANTLR start "ruleNumberValue"
    // InternalElkGraph.g:458:1: ruleNumberValue : ( ( rule__NumberValue__Alternatives ) ) ;
    public final void ruleNumberValue() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:462:2: ( ( ( rule__NumberValue__Alternatives ) ) )
            // InternalElkGraph.g:463:2: ( ( rule__NumberValue__Alternatives ) )
            {
            // InternalElkGraph.g:463:2: ( ( rule__NumberValue__Alternatives ) )
            // InternalElkGraph.g:464:3: ( rule__NumberValue__Alternatives )
            {
             before(grammarAccess.getNumberValueAccess().getAlternatives()); 
            // InternalElkGraph.g:465:3: ( rule__NumberValue__Alternatives )
            // InternalElkGraph.g:465:4: rule__NumberValue__Alternatives
            {
            pushFollow(FOLLOW_2);
            rule__NumberValue__Alternatives();

            state._fsp--;


            }

             after(grammarAccess.getNumberValueAccess().getAlternatives()); 

            }


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "ruleNumberValue"


    // $ANTLR start "entryRuleBooleanValue"
    // InternalElkGraph.g:474:1: entryRuleBooleanValue : ruleBooleanValue EOF ;
    public final void entryRuleBooleanValue() throws RecognitionException {
        try {
            // InternalElkGraph.g:475:1: ( ruleBooleanValue EOF )
            // InternalElkGraph.g:476:1: ruleBooleanValue EOF
            {
             before(grammarAccess.getBooleanValueRule()); 
            pushFollow(FOLLOW_1);
            ruleBooleanValue();

            state._fsp--;

             after(grammarAccess.getBooleanValueRule()); 
            match(input,EOF,FOLLOW_2); 

            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {
        }
        return ;
    }
    // $ANTLR end "entryRuleBooleanValue"


    // $ANTLR start "ruleBooleanValue"
    // InternalElkGraph.g:483:1: ruleBooleanValue : ( ( rule__BooleanValue__Alternatives ) ) ;
    public final void ruleBooleanValue() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:487:2: ( ( ( rule__BooleanValue__Alternatives ) ) )
            // InternalElkGraph.g:488:2: ( ( rule__BooleanValue__Alternatives ) )
            {
            // InternalElkGraph.g:488:2: ( ( rule__BooleanValue__Alternatives ) )
            // InternalElkGraph.g:489:3: ( rule__BooleanValue__Alternatives )
            {
             before(grammarAccess.getBooleanValueAccess().getAlternatives()); 
            // InternalElkGraph.g:490:3: ( rule__BooleanValue__Alternatives )
            // InternalElkGraph.g:490:4: rule__BooleanValue__Alternatives
            {
            pushFollow(FOLLOW_2);
            rule__BooleanValue__Alternatives();

            state._fsp--;


            }

             after(grammarAccess.getBooleanValueAccess().getAlternatives()); 

            }


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "ruleBooleanValue"


    // $ANTLR start "rule__RootNode__Alternatives_4"
    // InternalElkGraph.g:498:1: rule__RootNode__Alternatives_4 : ( ( ( rule__RootNode__LabelsAssignment_4_0 ) ) | ( ( rule__RootNode__PortsAssignment_4_1 ) ) | ( ( rule__RootNode__ChildrenAssignment_4_2 ) ) | ( ( rule__RootNode__ContainedEdgesAssignment_4_3 ) ) );
    public final void rule__RootNode__Alternatives_4() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:502:1: ( ( ( rule__RootNode__LabelsAssignment_4_0 ) ) | ( ( rule__RootNode__PortsAssignment_4_1 ) ) | ( ( rule__RootNode__ChildrenAssignment_4_2 ) ) | ( ( rule__RootNode__ContainedEdgesAssignment_4_3 ) ) )
            int alt1=4;
            switch ( input.LA(1) ) {
            case 20:
                {
                alt1=1;
                }
                break;
            case 22:
                {
                alt1=2;
                }
                break;
            case 17:
                {
                alt1=3;
                }
                break;
            case 29:
                {
                alt1=4;
                }
                break;
            default:
                NoViableAltException nvae =
                    new NoViableAltException("", 1, 0, input);

                throw nvae;
            }

            switch (alt1) {
                case 1 :
                    // InternalElkGraph.g:503:2: ( ( rule__RootNode__LabelsAssignment_4_0 ) )
                    {
                    // InternalElkGraph.g:503:2: ( ( rule__RootNode__LabelsAssignment_4_0 ) )
                    // InternalElkGraph.g:504:3: ( rule__RootNode__LabelsAssignment_4_0 )
                    {
                     before(grammarAccess.getRootNodeAccess().getLabelsAssignment_4_0()); 
                    // InternalElkGraph.g:505:3: ( rule__RootNode__LabelsAssignment_4_0 )
                    // InternalElkGraph.g:505:4: rule__RootNode__LabelsAssignment_4_0
                    {
                    pushFollow(FOLLOW_2);
                    rule__RootNode__LabelsAssignment_4_0();

                    state._fsp--;


                    }

                     after(grammarAccess.getRootNodeAccess().getLabelsAssignment_4_0()); 

                    }


                    }
                    break;
                case 2 :
                    // InternalElkGraph.g:509:2: ( ( rule__RootNode__PortsAssignment_4_1 ) )
                    {
                    // InternalElkGraph.g:509:2: ( ( rule__RootNode__PortsAssignment_4_1 ) )
                    // InternalElkGraph.g:510:3: ( rule__RootNode__PortsAssignment_4_1 )
                    {
                     before(grammarAccess.getRootNodeAccess().getPortsAssignment_4_1()); 
                    // InternalElkGraph.g:511:3: ( rule__RootNode__PortsAssignment_4_1 )
                    // InternalElkGraph.g:511:4: rule__RootNode__PortsAssignment_4_1
                    {
                    pushFollow(FOLLOW_2);
                    rule__RootNode__PortsAssignment_4_1();

                    state._fsp--;


                    }

                     after(grammarAccess.getRootNodeAccess().getPortsAssignment_4_1()); 

                    }


                    }
                    break;
                case 3 :
                    // InternalElkGraph.g:515:2: ( ( rule__RootNode__ChildrenAssignment_4_2 ) )
                    {
                    // InternalElkGraph.g:515:2: ( ( rule__RootNode__ChildrenAssignment_4_2 ) )
                    // InternalElkGraph.g:516:3: ( rule__RootNode__ChildrenAssignment_4_2 )
                    {
                     before(grammarAccess.getRootNodeAccess().getChildrenAssignment_4_2()); 
                    // InternalElkGraph.g:517:3: ( rule__RootNode__ChildrenAssignment_4_2 )
                    // InternalElkGraph.g:517:4: rule__RootNode__ChildrenAssignment_4_2
                    {
                    pushFollow(FOLLOW_2);
                    rule__RootNode__ChildrenAssignment_4_2();

                    state._fsp--;


                    }

                     after(grammarAccess.getRootNodeAccess().getChildrenAssignment_4_2()); 

                    }


                    }
                    break;
                case 4 :
                    // InternalElkGraph.g:521:2: ( ( rule__RootNode__ContainedEdgesAssignment_4_3 ) )
                    {
                    // InternalElkGraph.g:521:2: ( ( rule__RootNode__ContainedEdgesAssignment_4_3 ) )
                    // InternalElkGraph.g:522:3: ( rule__RootNode__ContainedEdgesAssignment_4_3 )
                    {
                     before(grammarAccess.getRootNodeAccess().getContainedEdgesAssignment_4_3()); 
                    // InternalElkGraph.g:523:3: ( rule__RootNode__ContainedEdgesAssignment_4_3 )
                    // InternalElkGraph.g:523:4: rule__RootNode__ContainedEdgesAssignment_4_3
                    {
                    pushFollow(FOLLOW_2);
                    rule__RootNode__ContainedEdgesAssignment_4_3();

                    state._fsp--;


                    }

                     after(grammarAccess.getRootNodeAccess().getContainedEdgesAssignment_4_3()); 

                    }


                    }
                    break;

            }
        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__RootNode__Alternatives_4"


    // $ANTLR start "rule__ElkNode__Alternatives_2_3"
    // InternalElkGraph.g:531:1: rule__ElkNode__Alternatives_2_3 : ( ( ( rule__ElkNode__LabelsAssignment_2_3_0 ) ) | ( ( rule__ElkNode__PortsAssignment_2_3_1 ) ) | ( ( rule__ElkNode__ChildrenAssignment_2_3_2 ) ) | ( ( rule__ElkNode__ContainedEdgesAssignment_2_3_3 ) ) );
    public final void rule__ElkNode__Alternatives_2_3() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:535:1: ( ( ( rule__ElkNode__LabelsAssignment_2_3_0 ) ) | ( ( rule__ElkNode__PortsAssignment_2_3_1 ) ) | ( ( rule__ElkNode__ChildrenAssignment_2_3_2 ) ) | ( ( rule__ElkNode__ContainedEdgesAssignment_2_3_3 ) ) )
            int alt2=4;
            switch ( input.LA(1) ) {
            case 20:
                {
                alt2=1;
                }
                break;
            case 22:
                {
                alt2=2;
                }
                break;
            case 17:
                {
                alt2=3;
                }
                break;
            case 29:
                {
                alt2=4;
                }
                break;
            default:
                NoViableAltException nvae =
                    new NoViableAltException("", 2, 0, input);

                throw nvae;
            }

            switch (alt2) {
                case 1 :
                    // InternalElkGraph.g:536:2: ( ( rule__ElkNode__LabelsAssignment_2_3_0 ) )
                    {
                    // InternalElkGraph.g:536:2: ( ( rule__ElkNode__LabelsAssignment_2_3_0 ) )
                    // InternalElkGraph.g:537:3: ( rule__ElkNode__LabelsAssignment_2_3_0 )
                    {
                     before(grammarAccess.getElkNodeAccess().getLabelsAssignment_2_3_0()); 
                    // InternalElkGraph.g:538:3: ( rule__ElkNode__LabelsAssignment_2_3_0 )
                    // InternalElkGraph.g:538:4: rule__ElkNode__LabelsAssignment_2_3_0
                    {
                    pushFollow(FOLLOW_2);
                    rule__ElkNode__LabelsAssignment_2_3_0();

                    state._fsp--;


                    }

                     after(grammarAccess.getElkNodeAccess().getLabelsAssignment_2_3_0()); 

                    }


                    }
                    break;
                case 2 :
                    // InternalElkGraph.g:542:2: ( ( rule__ElkNode__PortsAssignment_2_3_1 ) )
                    {
                    // InternalElkGraph.g:542:2: ( ( rule__ElkNode__PortsAssignment_2_3_1 ) )
                    // InternalElkGraph.g:543:3: ( rule__ElkNode__PortsAssignment_2_3_1 )
                    {
                     before(grammarAccess.getElkNodeAccess().getPortsAssignment_2_3_1()); 
                    // InternalElkGraph.g:544:3: ( rule__ElkNode__PortsAssignment_2_3_1 )
                    // InternalElkGraph.g:544:4: rule__ElkNode__PortsAssignment_2_3_1
                    {
                    pushFollow(FOLLOW_2);
                    rule__ElkNode__PortsAssignment_2_3_1();

                    state._fsp--;


                    }

                     after(grammarAccess.getElkNodeAccess().getPortsAssignment_2_3_1()); 

                    }


                    }
                    break;
                case 3 :
                    // InternalElkGraph.g:548:2: ( ( rule__ElkNode__ChildrenAssignment_2_3_2 ) )
                    {
                    // InternalElkGraph.g:548:2: ( ( rule__ElkNode__ChildrenAssignment_2_3_2 ) )
                    // InternalElkGraph.g:549:3: ( rule__ElkNode__ChildrenAssignment_2_3_2 )
                    {
                     before(grammarAccess.getElkNodeAccess().getChildrenAssignment_2_3_2()); 
                    // InternalElkGraph.g:550:3: ( rule__ElkNode__ChildrenAssignment_2_3_2 )
                    // InternalElkGraph.g:550:4: rule__ElkNode__ChildrenAssignment_2_3_2
                    {
                    pushFollow(FOLLOW_2);
                    rule__ElkNode__ChildrenAssignment_2_3_2();

                    state._fsp--;


                    }

                     after(grammarAccess.getElkNodeAccess().getChildrenAssignment_2_3_2()); 

                    }


                    }
                    break;
                case 4 :
                    // InternalElkGraph.g:554:2: ( ( rule__ElkNode__ContainedEdgesAssignment_2_3_3 ) )
                    {
                    // InternalElkGraph.g:554:2: ( ( rule__ElkNode__ContainedEdgesAssignment_2_3_3 ) )
                    // InternalElkGraph.g:555:3: ( rule__ElkNode__ContainedEdgesAssignment_2_3_3 )
                    {
                     before(grammarAccess.getElkNodeAccess().getContainedEdgesAssignment_2_3_3()); 
                    // InternalElkGraph.g:556:3: ( rule__ElkNode__ContainedEdgesAssignment_2_3_3 )
                    // InternalElkGraph.g:556:4: rule__ElkNode__ContainedEdgesAssignment_2_3_3
                    {
                    pushFollow(FOLLOW_2);
                    rule__ElkNode__ContainedEdgesAssignment_2_3_3();

                    state._fsp--;


                    }

                     after(grammarAccess.getElkNodeAccess().getContainedEdgesAssignment_2_3_3()); 

                    }


                    }
                    break;

            }
        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ElkNode__Alternatives_2_3"


    // $ANTLR start "rule__EdgeLayout__Alternatives_2"
    // InternalElkGraph.g:564:1: rule__EdgeLayout__Alternatives_2 : ( ( ( rule__EdgeLayout__SectionsAssignment_2_0 ) ) | ( ( ( rule__EdgeLayout__SectionsAssignment_2_1 ) ) ( ( rule__EdgeLayout__SectionsAssignment_2_1 )* ) ) );
    public final void rule__EdgeLayout__Alternatives_2() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:568:1: ( ( ( rule__EdgeLayout__SectionsAssignment_2_0 ) ) | ( ( ( rule__EdgeLayout__SectionsAssignment_2_1 ) ) ( ( rule__EdgeLayout__SectionsAssignment_2_1 )* ) ) )
            int alt4=2;
            int LA4_0 = input.LA(1);

            if ( (LA4_0==EOF||LA4_0==RULE_ID||LA4_0==25||(LA4_0>=31 && LA4_0<=35)) ) {
                alt4=1;
            }
            else if ( (LA4_0==37) ) {
                alt4=2;
            }
            else {
                NoViableAltException nvae =
                    new NoViableAltException("", 4, 0, input);

                throw nvae;
            }
            switch (alt4) {
                case 1 :
                    // InternalElkGraph.g:569:2: ( ( rule__EdgeLayout__SectionsAssignment_2_0 ) )
                    {
                    // InternalElkGraph.g:569:2: ( ( rule__EdgeLayout__SectionsAssignment_2_0 ) )
                    // InternalElkGraph.g:570:3: ( rule__EdgeLayout__SectionsAssignment_2_0 )
                    {
                     before(grammarAccess.getEdgeLayoutAccess().getSectionsAssignment_2_0()); 
                    // InternalElkGraph.g:571:3: ( rule__EdgeLayout__SectionsAssignment_2_0 )
                    // InternalElkGraph.g:571:4: rule__EdgeLayout__SectionsAssignment_2_0
                    {
                    pushFollow(FOLLOW_2);
                    rule__EdgeLayout__SectionsAssignment_2_0();

                    state._fsp--;


                    }

                     after(grammarAccess.getEdgeLayoutAccess().getSectionsAssignment_2_0()); 

                    }


                    }
                    break;
                case 2 :
                    // InternalElkGraph.g:575:2: ( ( ( rule__EdgeLayout__SectionsAssignment_2_1 ) ) ( ( rule__EdgeLayout__SectionsAssignment_2_1 )* ) )
                    {
                    // InternalElkGraph.g:575:2: ( ( ( rule__EdgeLayout__SectionsAssignment_2_1 ) ) ( ( rule__EdgeLayout__SectionsAssignment_2_1 )* ) )
                    // InternalElkGraph.g:576:3: ( ( rule__EdgeLayout__SectionsAssignment_2_1 ) ) ( ( rule__EdgeLayout__SectionsAssignment_2_1 )* )
                    {
                    // InternalElkGraph.g:576:3: ( ( rule__EdgeLayout__SectionsAssignment_2_1 ) )
                    // InternalElkGraph.g:577:4: ( rule__EdgeLayout__SectionsAssignment_2_1 )
                    {
                     before(grammarAccess.getEdgeLayoutAccess().getSectionsAssignment_2_1()); 
                    // InternalElkGraph.g:578:4: ( rule__EdgeLayout__SectionsAssignment_2_1 )
                    // InternalElkGraph.g:578:5: rule__EdgeLayout__SectionsAssignment_2_1
                    {
                    pushFollow(FOLLOW_3);
                    rule__EdgeLayout__SectionsAssignment_2_1();

                    state._fsp--;


                    }

                     after(grammarAccess.getEdgeLayoutAccess().getSectionsAssignment_2_1()); 

                    }

                    // InternalElkGraph.g:581:3: ( ( rule__EdgeLayout__SectionsAssignment_2_1 )* )
                    // InternalElkGraph.g:582:4: ( rule__EdgeLayout__SectionsAssignment_2_1 )*
                    {
                     before(grammarAccess.getEdgeLayoutAccess().getSectionsAssignment_2_1()); 
                    // InternalElkGraph.g:583:4: ( rule__EdgeLayout__SectionsAssignment_2_1 )*
                    loop3:
                    do {
                        int alt3=2;
                        int LA3_0 = input.LA(1);

                        if ( (LA3_0==37) ) {
                            alt3=1;
                        }


                        switch (alt3) {
                    	case 1 :
                    	    // InternalElkGraph.g:583:5: rule__EdgeLayout__SectionsAssignment_2_1
                    	    {
                    	    pushFollow(FOLLOW_3);
                    	    rule__EdgeLayout__SectionsAssignment_2_1();

                    	    state._fsp--;


                    	    }
                    	    break;

                    	default :
                    	    break loop3;
                        }
                    } while (true);

                     after(grammarAccess.getEdgeLayoutAccess().getSectionsAssignment_2_1()); 

                    }


                    }


                    }
                    break;

            }
        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__EdgeLayout__Alternatives_2"


    // $ANTLR start "rule__Number__Alternatives"
    // InternalElkGraph.g:592:1: rule__Number__Alternatives : ( ( RULE_SIGNED_INT ) | ( RULE_FLOAT ) );
    public final void rule__Number__Alternatives() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:596:1: ( ( RULE_SIGNED_INT ) | ( RULE_FLOAT ) )
            int alt5=2;
            int LA5_0 = input.LA(1);

            if ( (LA5_0==RULE_SIGNED_INT) ) {
                alt5=1;
            }
            else if ( (LA5_0==RULE_FLOAT) ) {
                alt5=2;
            }
            else {
                NoViableAltException nvae =
                    new NoViableAltException("", 5, 0, input);

                throw nvae;
            }
            switch (alt5) {
                case 1 :
                    // InternalElkGraph.g:597:2: ( RULE_SIGNED_INT )
                    {
                    // InternalElkGraph.g:597:2: ( RULE_SIGNED_INT )
                    // InternalElkGraph.g:598:3: RULE_SIGNED_INT
                    {
                     before(grammarAccess.getNumberAccess().getSIGNED_INTTerminalRuleCall_0()); 
                    match(input,RULE_SIGNED_INT,FOLLOW_2); 
                     after(grammarAccess.getNumberAccess().getSIGNED_INTTerminalRuleCall_0()); 

                    }


                    }
                    break;
                case 2 :
                    // InternalElkGraph.g:603:2: ( RULE_FLOAT )
                    {
                    // InternalElkGraph.g:603:2: ( RULE_FLOAT )
                    // InternalElkGraph.g:604:3: RULE_FLOAT
                    {
                     before(grammarAccess.getNumberAccess().getFLOATTerminalRuleCall_1()); 
                    match(input,RULE_FLOAT,FOLLOW_2); 
                     after(grammarAccess.getNumberAccess().getFLOATTerminalRuleCall_1()); 

                    }


                    }
                    break;

            }
        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__Number__Alternatives"


    // $ANTLR start "rule__Property__Alternatives_2"
    // InternalElkGraph.g:613:1: rule__Property__Alternatives_2 : ( ( ( rule__Property__ValueAssignment_2_0 ) ) | ( ( rule__Property__ValueAssignment_2_1 ) ) | ( ( rule__Property__ValueAssignment_2_2 ) ) | ( ( rule__Property__ValueAssignment_2_3 ) ) | ( 'null' ) );
    public final void rule__Property__Alternatives_2() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:617:1: ( ( ( rule__Property__ValueAssignment_2_0 ) ) | ( ( rule__Property__ValueAssignment_2_1 ) ) | ( ( rule__Property__ValueAssignment_2_2 ) ) | ( ( rule__Property__ValueAssignment_2_3 ) ) | ( 'null' ) )
            int alt6=5;
            switch ( input.LA(1) ) {
            case RULE_STRING:
                {
                alt6=1;
                }
                break;
            case RULE_ID:
                {
                alt6=2;
                }
                break;
            case RULE_SIGNED_INT:
            case RULE_FLOAT:
                {
                alt6=3;
                }
                break;
            case 14:
            case 15:
                {
                alt6=4;
                }
                break;
            case 13:
                {
                alt6=5;
                }
                break;
            default:
                NoViableAltException nvae =
                    new NoViableAltException("", 6, 0, input);

                throw nvae;
            }

            switch (alt6) {
                case 1 :
                    // InternalElkGraph.g:618:2: ( ( rule__Property__ValueAssignment_2_0 ) )
                    {
                    // InternalElkGraph.g:618:2: ( ( rule__Property__ValueAssignment_2_0 ) )
                    // InternalElkGraph.g:619:3: ( rule__Property__ValueAssignment_2_0 )
                    {
                     before(grammarAccess.getPropertyAccess().getValueAssignment_2_0()); 
                    // InternalElkGraph.g:620:3: ( rule__Property__ValueAssignment_2_0 )
                    // InternalElkGraph.g:620:4: rule__Property__ValueAssignment_2_0
                    {
                    pushFollow(FOLLOW_2);
                    rule__Property__ValueAssignment_2_0();

                    state._fsp--;


                    }

                     after(grammarAccess.getPropertyAccess().getValueAssignment_2_0()); 

                    }


                    }
                    break;
                case 2 :
                    // InternalElkGraph.g:624:2: ( ( rule__Property__ValueAssignment_2_1 ) )
                    {
                    // InternalElkGraph.g:624:2: ( ( rule__Property__ValueAssignment_2_1 ) )
                    // InternalElkGraph.g:625:3: ( rule__Property__ValueAssignment_2_1 )
                    {
                     before(grammarAccess.getPropertyAccess().getValueAssignment_2_1()); 
                    // InternalElkGraph.g:626:3: ( rule__Property__ValueAssignment_2_1 )
                    // InternalElkGraph.g:626:4: rule__Property__ValueAssignment_2_1
                    {
                    pushFollow(FOLLOW_2);
                    rule__Property__ValueAssignment_2_1();

                    state._fsp--;


                    }

                     after(grammarAccess.getPropertyAccess().getValueAssignment_2_1()); 

                    }


                    }
                    break;
                case 3 :
                    // InternalElkGraph.g:630:2: ( ( rule__Property__ValueAssignment_2_2 ) )
                    {
                    // InternalElkGraph.g:630:2: ( ( rule__Property__ValueAssignment_2_2 ) )
                    // InternalElkGraph.g:631:3: ( rule__Property__ValueAssignment_2_2 )
                    {
                     before(grammarAccess.getPropertyAccess().getValueAssignment_2_2()); 
                    // InternalElkGraph.g:632:3: ( rule__Property__ValueAssignment_2_2 )
                    // InternalElkGraph.g:632:4: rule__Property__ValueAssignment_2_2
                    {
                    pushFollow(FOLLOW_2);
                    rule__Property__ValueAssignment_2_2();

                    state._fsp--;


                    }

                     after(grammarAccess.getPropertyAccess().getValueAssignment_2_2()); 

                    }


                    }
                    break;
                case 4 :
                    // InternalElkGraph.g:636:2: ( ( rule__Property__ValueAssignment_2_3 ) )
                    {
                    // InternalElkGraph.g:636:2: ( ( rule__Property__ValueAssignment_2_3 ) )
                    // InternalElkGraph.g:637:3: ( rule__Property__ValueAssignment_2_3 )
                    {
                     before(grammarAccess.getPropertyAccess().getValueAssignment_2_3()); 
                    // InternalElkGraph.g:638:3: ( rule__Property__ValueAssignment_2_3 )
                    // InternalElkGraph.g:638:4: rule__Property__ValueAssignment_2_3
                    {
                    pushFollow(FOLLOW_2);
                    rule__Property__ValueAssignment_2_3();

                    state._fsp--;


                    }

                     after(grammarAccess.getPropertyAccess().getValueAssignment_2_3()); 

                    }


                    }
                    break;
                case 5 :
                    // InternalElkGraph.g:642:2: ( 'null' )
                    {
                    // InternalElkGraph.g:642:2: ( 'null' )
                    // InternalElkGraph.g:643:3: 'null'
                    {
                     before(grammarAccess.getPropertyAccess().getNullKeyword_2_4()); 
                    match(input,13,FOLLOW_2); 
                     after(grammarAccess.getPropertyAccess().getNullKeyword_2_4()); 

                    }


                    }
                    break;

            }
        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__Property__Alternatives_2"


    // $ANTLR start "rule__NumberValue__Alternatives"
    // InternalElkGraph.g:652:1: rule__NumberValue__Alternatives : ( ( RULE_SIGNED_INT ) | ( RULE_FLOAT ) );
    public final void rule__NumberValue__Alternatives() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:656:1: ( ( RULE_SIGNED_INT ) | ( RULE_FLOAT ) )
            int alt7=2;
            int LA7_0 = input.LA(1);

            if ( (LA7_0==RULE_SIGNED_INT) ) {
                alt7=1;
            }
            else if ( (LA7_0==RULE_FLOAT) ) {
                alt7=2;
            }
            else {
                NoViableAltException nvae =
                    new NoViableAltException("", 7, 0, input);

                throw nvae;
            }
            switch (alt7) {
                case 1 :
                    // InternalElkGraph.g:657:2: ( RULE_SIGNED_INT )
                    {
                    // InternalElkGraph.g:657:2: ( RULE_SIGNED_INT )
                    // InternalElkGraph.g:658:3: RULE_SIGNED_INT
                    {
                     before(grammarAccess.getNumberValueAccess().getSIGNED_INTTerminalRuleCall_0()); 
                    match(input,RULE_SIGNED_INT,FOLLOW_2); 
                     after(grammarAccess.getNumberValueAccess().getSIGNED_INTTerminalRuleCall_0()); 

                    }


                    }
                    break;
                case 2 :
                    // InternalElkGraph.g:663:2: ( RULE_FLOAT )
                    {
                    // InternalElkGraph.g:663:2: ( RULE_FLOAT )
                    // InternalElkGraph.g:664:3: RULE_FLOAT
                    {
                     before(grammarAccess.getNumberValueAccess().getFLOATTerminalRuleCall_1()); 
                    match(input,RULE_FLOAT,FOLLOW_2); 
                     after(grammarAccess.getNumberValueAccess().getFLOATTerminalRuleCall_1()); 

                    }


                    }
                    break;

            }
        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__NumberValue__Alternatives"


    // $ANTLR start "rule__BooleanValue__Alternatives"
    // InternalElkGraph.g:673:1: rule__BooleanValue__Alternatives : ( ( 'true' ) | ( 'false' ) );
    public final void rule__BooleanValue__Alternatives() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:677:1: ( ( 'true' ) | ( 'false' ) )
            int alt8=2;
            int LA8_0 = input.LA(1);

            if ( (LA8_0==14) ) {
                alt8=1;
            }
            else if ( (LA8_0==15) ) {
                alt8=2;
            }
            else {
                NoViableAltException nvae =
                    new NoViableAltException("", 8, 0, input);

                throw nvae;
            }
            switch (alt8) {
                case 1 :
                    // InternalElkGraph.g:678:2: ( 'true' )
                    {
                    // InternalElkGraph.g:678:2: ( 'true' )
                    // InternalElkGraph.g:679:3: 'true'
                    {
                     before(grammarAccess.getBooleanValueAccess().getTrueKeyword_0()); 
                    match(input,14,FOLLOW_2); 
                     after(grammarAccess.getBooleanValueAccess().getTrueKeyword_0()); 

                    }


                    }
                    break;
                case 2 :
                    // InternalElkGraph.g:684:2: ( 'false' )
                    {
                    // InternalElkGraph.g:684:2: ( 'false' )
                    // InternalElkGraph.g:685:3: 'false'
                    {
                     before(grammarAccess.getBooleanValueAccess().getFalseKeyword_1()); 
                    match(input,15,FOLLOW_2); 
                     after(grammarAccess.getBooleanValueAccess().getFalseKeyword_1()); 

                    }


                    }
                    break;

            }
        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__BooleanValue__Alternatives"


    // $ANTLR start "rule__RootNode__Group__0"
    // InternalElkGraph.g:694:1: rule__RootNode__Group__0 : rule__RootNode__Group__0__Impl rule__RootNode__Group__1 ;
    public final void rule__RootNode__Group__0() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:698:1: ( rule__RootNode__Group__0__Impl rule__RootNode__Group__1 )
            // InternalElkGraph.g:699:2: rule__RootNode__Group__0__Impl rule__RootNode__Group__1
            {
            pushFollow(FOLLOW_4);
            rule__RootNode__Group__0__Impl();

            state._fsp--;

            pushFollow(FOLLOW_2);
            rule__RootNode__Group__1();

            state._fsp--;


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__RootNode__Group__0"


    // $ANTLR start "rule__RootNode__Group__0__Impl"
    // InternalElkGraph.g:706:1: rule__RootNode__Group__0__Impl : ( () ) ;
    public final void rule__RootNode__Group__0__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:710:1: ( ( () ) )
            // InternalElkGraph.g:711:1: ( () )
            {
            // InternalElkGraph.g:711:1: ( () )
            // InternalElkGraph.g:712:2: ()
            {
             before(grammarAccess.getRootNodeAccess().getElkNodeAction_0()); 
            // InternalElkGraph.g:713:2: ()
            // InternalElkGraph.g:713:3: 
            {
            }

             after(grammarAccess.getRootNodeAccess().getElkNodeAction_0()); 

            }


            }

        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__RootNode__Group__0__Impl"


    // $ANTLR start "rule__RootNode__Group__1"
    // InternalElkGraph.g:721:1: rule__RootNode__Group__1 : rule__RootNode__Group__1__Impl rule__RootNode__Group__2 ;
    public final void rule__RootNode__Group__1() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:725:1: ( rule__RootNode__Group__1__Impl rule__RootNode__Group__2 )
            // InternalElkGraph.g:726:2: rule__RootNode__Group__1__Impl rule__RootNode__Group__2
            {
            pushFollow(FOLLOW_4);
            rule__RootNode__Group__1__Impl();

            state._fsp--;

            pushFollow(FOLLOW_2);
            rule__RootNode__Group__2();

            state._fsp--;


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__RootNode__Group__1"


    // $ANTLR start "rule__RootNode__Group__1__Impl"
    // InternalElkGraph.g:733:1: rule__RootNode__Group__1__Impl : ( ( rule__RootNode__Group_1__0 )? ) ;
    public final void rule__RootNode__Group__1__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:737:1: ( ( ( rule__RootNode__Group_1__0 )? ) )
            // InternalElkGraph.g:738:1: ( ( rule__RootNode__Group_1__0 )? )
            {
            // InternalElkGraph.g:738:1: ( ( rule__RootNode__Group_1__0 )? )
            // InternalElkGraph.g:739:2: ( rule__RootNode__Group_1__0 )?
            {
             before(grammarAccess.getRootNodeAccess().getGroup_1()); 
            // InternalElkGraph.g:740:2: ( rule__RootNode__Group_1__0 )?
            int alt9=2;
            int LA9_0 = input.LA(1);

            if ( (LA9_0==16) ) {
                alt9=1;
            }
            switch (alt9) {
                case 1 :
                    // InternalElkGraph.g:740:3: rule__RootNode__Group_1__0
                    {
                    pushFollow(FOLLOW_2);
                    rule__RootNode__Group_1__0();

                    state._fsp--;


                    }
                    break;

            }

             after(grammarAccess.getRootNodeAccess().getGroup_1()); 

            }


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__RootNode__Group__1__Impl"


    // $ANTLR start "rule__RootNode__Group__2"
    // InternalElkGraph.g:748:1: rule__RootNode__Group__2 : rule__RootNode__Group__2__Impl rule__RootNode__Group__3 ;
    public final void rule__RootNode__Group__2() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:752:1: ( rule__RootNode__Group__2__Impl rule__RootNode__Group__3 )
            // InternalElkGraph.g:753:2: rule__RootNode__Group__2__Impl rule__RootNode__Group__3
            {
            pushFollow(FOLLOW_4);
            rule__RootNode__Group__2__Impl();

            state._fsp--;

            pushFollow(FOLLOW_2);
            rule__RootNode__Group__3();

            state._fsp--;


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__RootNode__Group__2"


    // $ANTLR start "rule__RootNode__Group__2__Impl"
    // InternalElkGraph.g:760:1: rule__RootNode__Group__2__Impl : ( ( ruleShapeLayout )? ) ;
    public final void rule__RootNode__Group__2__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:764:1: ( ( ( ruleShapeLayout )? ) )
            // InternalElkGraph.g:765:1: ( ( ruleShapeLayout )? )
            {
            // InternalElkGraph.g:765:1: ( ( ruleShapeLayout )? )
            // InternalElkGraph.g:766:2: ( ruleShapeLayout )?
            {
             before(grammarAccess.getRootNodeAccess().getShapeLayoutParserRuleCall_2()); 
            // InternalElkGraph.g:767:2: ( ruleShapeLayout )?
            int alt10=2;
            int LA10_0 = input.LA(1);

            if ( (LA10_0==23) ) {
                alt10=1;
            }
            switch (alt10) {
                case 1 :
                    // InternalElkGraph.g:767:3: ruleShapeLayout
                    {
                    pushFollow(FOLLOW_2);
                    ruleShapeLayout();

                    state._fsp--;


                    }
                    break;

            }

             after(grammarAccess.getRootNodeAccess().getShapeLayoutParserRuleCall_2()); 

            }


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__RootNode__Group__2__Impl"


    // $ANTLR start "rule__RootNode__Group__3"
    // InternalElkGraph.g:775:1: rule__RootNode__Group__3 : rule__RootNode__Group__3__Impl rule__RootNode__Group__4 ;
    public final void rule__RootNode__Group__3() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:779:1: ( rule__RootNode__Group__3__Impl rule__RootNode__Group__4 )
            // InternalElkGraph.g:780:2: rule__RootNode__Group__3__Impl rule__RootNode__Group__4
            {
            pushFollow(FOLLOW_4);
            rule__RootNode__Group__3__Impl();

            state._fsp--;

            pushFollow(FOLLOW_2);
            rule__RootNode__Group__4();

            state._fsp--;


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__RootNode__Group__3"


    // $ANTLR start "rule__RootNode__Group__3__Impl"
    // InternalElkGraph.g:787:1: rule__RootNode__Group__3__Impl : ( ( rule__RootNode__PropertiesAssignment_3 )* ) ;
    public final void rule__RootNode__Group__3__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:791:1: ( ( ( rule__RootNode__PropertiesAssignment_3 )* ) )
            // InternalElkGraph.g:792:1: ( ( rule__RootNode__PropertiesAssignment_3 )* )
            {
            // InternalElkGraph.g:792:1: ( ( rule__RootNode__PropertiesAssignment_3 )* )
            // InternalElkGraph.g:793:2: ( rule__RootNode__PropertiesAssignment_3 )*
            {
             before(grammarAccess.getRootNodeAccess().getPropertiesAssignment_3()); 
            // InternalElkGraph.g:794:2: ( rule__RootNode__PropertiesAssignment_3 )*
            loop11:
            do {
                int alt11=2;
                int LA11_0 = input.LA(1);

                if ( (LA11_0==RULE_ID) ) {
                    alt11=1;
                }


                switch (alt11) {
            	case 1 :
            	    // InternalElkGraph.g:794:3: rule__RootNode__PropertiesAssignment_3
            	    {
            	    pushFollow(FOLLOW_5);
            	    rule__RootNode__PropertiesAssignment_3();

            	    state._fsp--;


            	    }
            	    break;

            	default :
            	    break loop11;
                }
            } while (true);

             after(grammarAccess.getRootNodeAccess().getPropertiesAssignment_3()); 

            }


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__RootNode__Group__3__Impl"


    // $ANTLR start "rule__RootNode__Group__4"
    // InternalElkGraph.g:802:1: rule__RootNode__Group__4 : rule__RootNode__Group__4__Impl ;
    public final void rule__RootNode__Group__4() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:806:1: ( rule__RootNode__Group__4__Impl )
            // InternalElkGraph.g:807:2: rule__RootNode__Group__4__Impl
            {
            pushFollow(FOLLOW_2);
            rule__RootNode__Group__4__Impl();

            state._fsp--;


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__RootNode__Group__4"


    // $ANTLR start "rule__RootNode__Group__4__Impl"
    // InternalElkGraph.g:813:1: rule__RootNode__Group__4__Impl : ( ( rule__RootNode__Alternatives_4 )* ) ;
    public final void rule__RootNode__Group__4__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:817:1: ( ( ( rule__RootNode__Alternatives_4 )* ) )
            // InternalElkGraph.g:818:1: ( ( rule__RootNode__Alternatives_4 )* )
            {
            // InternalElkGraph.g:818:1: ( ( rule__RootNode__Alternatives_4 )* )
            // InternalElkGraph.g:819:2: ( rule__RootNode__Alternatives_4 )*
            {
             before(grammarAccess.getRootNodeAccess().getAlternatives_4()); 
            // InternalElkGraph.g:820:2: ( rule__RootNode__Alternatives_4 )*
            loop12:
            do {
                int alt12=2;
                int LA12_0 = input.LA(1);

                if ( (LA12_0==17||LA12_0==20||LA12_0==22||LA12_0==29) ) {
                    alt12=1;
                }


                switch (alt12) {
            	case 1 :
            	    // InternalElkGraph.g:820:3: rule__RootNode__Alternatives_4
            	    {
            	    pushFollow(FOLLOW_6);
            	    rule__RootNode__Alternatives_4();

            	    state._fsp--;


            	    }
            	    break;

            	default :
            	    break loop12;
                }
            } while (true);

             after(grammarAccess.getRootNodeAccess().getAlternatives_4()); 

            }


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__RootNode__Group__4__Impl"


    // $ANTLR start "rule__RootNode__Group_1__0"
    // InternalElkGraph.g:829:1: rule__RootNode__Group_1__0 : rule__RootNode__Group_1__0__Impl rule__RootNode__Group_1__1 ;
    public final void rule__RootNode__Group_1__0() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:833:1: ( rule__RootNode__Group_1__0__Impl rule__RootNode__Group_1__1 )
            // InternalElkGraph.g:834:2: rule__RootNode__Group_1__0__Impl rule__RootNode__Group_1__1
            {
            pushFollow(FOLLOW_7);
            rule__RootNode__Group_1__0__Impl();

            state._fsp--;

            pushFollow(FOLLOW_2);
            rule__RootNode__Group_1__1();

            state._fsp--;


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__RootNode__Group_1__0"


    // $ANTLR start "rule__RootNode__Group_1__0__Impl"
    // InternalElkGraph.g:841:1: rule__RootNode__Group_1__0__Impl : ( 'graph' ) ;
    public final void rule__RootNode__Group_1__0__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:845:1: ( ( 'graph' ) )
            // InternalElkGraph.g:846:1: ( 'graph' )
            {
            // InternalElkGraph.g:846:1: ( 'graph' )
            // InternalElkGraph.g:847:2: 'graph'
            {
             before(grammarAccess.getRootNodeAccess().getGraphKeyword_1_0()); 
            match(input,16,FOLLOW_2); 
             after(grammarAccess.getRootNodeAccess().getGraphKeyword_1_0()); 

            }


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__RootNode__Group_1__0__Impl"


    // $ANTLR start "rule__RootNode__Group_1__1"
    // InternalElkGraph.g:856:1: rule__RootNode__Group_1__1 : rule__RootNode__Group_1__1__Impl ;
    public final void rule__RootNode__Group_1__1() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:860:1: ( rule__RootNode__Group_1__1__Impl )
            // InternalElkGraph.g:861:2: rule__RootNode__Group_1__1__Impl
            {
            pushFollow(FOLLOW_2);
            rule__RootNode__Group_1__1__Impl();

            state._fsp--;


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__RootNode__Group_1__1"


    // $ANTLR start "rule__RootNode__Group_1__1__Impl"
    // InternalElkGraph.g:867:1: rule__RootNode__Group_1__1__Impl : ( ( rule__RootNode__IdentifierAssignment_1_1 ) ) ;
    public final void rule__RootNode__Group_1__1__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:871:1: ( ( ( rule__RootNode__IdentifierAssignment_1_1 ) ) )
            // InternalElkGraph.g:872:1: ( ( rule__RootNode__IdentifierAssignment_1_1 ) )
            {
            // InternalElkGraph.g:872:1: ( ( rule__RootNode__IdentifierAssignment_1_1 ) )
            // InternalElkGraph.g:873:2: ( rule__RootNode__IdentifierAssignment_1_1 )
            {
             before(grammarAccess.getRootNodeAccess().getIdentifierAssignment_1_1()); 
            // InternalElkGraph.g:874:2: ( rule__RootNode__IdentifierAssignment_1_1 )
            // InternalElkGraph.g:874:3: rule__RootNode__IdentifierAssignment_1_1
            {
            pushFollow(FOLLOW_2);
            rule__RootNode__IdentifierAssignment_1_1();

            state._fsp--;


            }

             after(grammarAccess.getRootNodeAccess().getIdentifierAssignment_1_1()); 

            }


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__RootNode__Group_1__1__Impl"


    // $ANTLR start "rule__ElkNode__Group__0"
    // InternalElkGraph.g:883:1: rule__ElkNode__Group__0 : rule__ElkNode__Group__0__Impl rule__ElkNode__Group__1 ;
    public final void rule__ElkNode__Group__0() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:887:1: ( rule__ElkNode__Group__0__Impl rule__ElkNode__Group__1 )
            // InternalElkGraph.g:888:2: rule__ElkNode__Group__0__Impl rule__ElkNode__Group__1
            {
            pushFollow(FOLLOW_7);
            rule__ElkNode__Group__0__Impl();

            state._fsp--;

            pushFollow(FOLLOW_2);
            rule__ElkNode__Group__1();

            state._fsp--;


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ElkNode__Group__0"


    // $ANTLR start "rule__ElkNode__Group__0__Impl"
    // InternalElkGraph.g:895:1: rule__ElkNode__Group__0__Impl : ( 'node' ) ;
    public final void rule__ElkNode__Group__0__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:899:1: ( ( 'node' ) )
            // InternalElkGraph.g:900:1: ( 'node' )
            {
            // InternalElkGraph.g:900:1: ( 'node' )
            // InternalElkGraph.g:901:2: 'node'
            {
             before(grammarAccess.getElkNodeAccess().getNodeKeyword_0()); 
            match(input,17,FOLLOW_2); 
             after(grammarAccess.getElkNodeAccess().getNodeKeyword_0()); 

            }


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ElkNode__Group__0__Impl"


    // $ANTLR start "rule__ElkNode__Group__1"
    // InternalElkGraph.g:910:1: rule__ElkNode__Group__1 : rule__ElkNode__Group__1__Impl rule__ElkNode__Group__2 ;
    public final void rule__ElkNode__Group__1() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:914:1: ( rule__ElkNode__Group__1__Impl rule__ElkNode__Group__2 )
            // InternalElkGraph.g:915:2: rule__ElkNode__Group__1__Impl rule__ElkNode__Group__2
            {
            pushFollow(FOLLOW_8);
            rule__ElkNode__Group__1__Impl();

            state._fsp--;

            pushFollow(FOLLOW_2);
            rule__ElkNode__Group__2();

            state._fsp--;


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ElkNode__Group__1"


    // $ANTLR start "rule__ElkNode__Group__1__Impl"
    // InternalElkGraph.g:922:1: rule__ElkNode__Group__1__Impl : ( ( rule__ElkNode__IdentifierAssignment_1 ) ) ;
    public final void rule__ElkNode__Group__1__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:926:1: ( ( ( rule__ElkNode__IdentifierAssignment_1 ) ) )
            // InternalElkGraph.g:927:1: ( ( rule__ElkNode__IdentifierAssignment_1 ) )
            {
            // InternalElkGraph.g:927:1: ( ( rule__ElkNode__IdentifierAssignment_1 ) )
            // InternalElkGraph.g:928:2: ( rule__ElkNode__IdentifierAssignment_1 )
            {
             before(grammarAccess.getElkNodeAccess().getIdentifierAssignment_1()); 
            // InternalElkGraph.g:929:2: ( rule__ElkNode__IdentifierAssignment_1 )
            // InternalElkGraph.g:929:3: rule__ElkNode__IdentifierAssignment_1
            {
            pushFollow(FOLLOW_2);
            rule__ElkNode__IdentifierAssignment_1();

            state._fsp--;


            }

             after(grammarAccess.getElkNodeAccess().getIdentifierAssignment_1()); 

            }


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ElkNode__Group__1__Impl"


    // $ANTLR start "rule__ElkNode__Group__2"
    // InternalElkGraph.g:937:1: rule__ElkNode__Group__2 : rule__ElkNode__Group__2__Impl ;
    public final void rule__ElkNode__Group__2() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:941:1: ( rule__ElkNode__Group__2__Impl )
            // InternalElkGraph.g:942:2: rule__ElkNode__Group__2__Impl
            {
            pushFollow(FOLLOW_2);
            rule__ElkNode__Group__2__Impl();

            state._fsp--;


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ElkNode__Group__2"


    // $ANTLR start "rule__ElkNode__Group__2__Impl"
    // InternalElkGraph.g:948:1: rule__ElkNode__Group__2__Impl : ( ( rule__ElkNode__Group_2__0 )? ) ;
    public final void rule__ElkNode__Group__2__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:952:1: ( ( ( rule__ElkNode__Group_2__0 )? ) )
            // InternalElkGraph.g:953:1: ( ( rule__ElkNode__Group_2__0 )? )
            {
            // InternalElkGraph.g:953:1: ( ( rule__ElkNode__Group_2__0 )? )
            // InternalElkGraph.g:954:2: ( rule__ElkNode__Group_2__0 )?
            {
             before(grammarAccess.getElkNodeAccess().getGroup_2()); 
            // InternalElkGraph.g:955:2: ( rule__ElkNode__Group_2__0 )?
            int alt13=2;
            int LA13_0 = input.LA(1);

            if ( (LA13_0==18) ) {
                alt13=1;
            }
            switch (alt13) {
                case 1 :
                    // InternalElkGraph.g:955:3: rule__ElkNode__Group_2__0
                    {
                    pushFollow(FOLLOW_2);
                    rule__ElkNode__Group_2__0();

                    state._fsp--;


                    }
                    break;

            }

             after(grammarAccess.getElkNodeAccess().getGroup_2()); 

            }


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ElkNode__Group__2__Impl"


    // $ANTLR start "rule__ElkNode__Group_2__0"
    // InternalElkGraph.g:964:1: rule__ElkNode__Group_2__0 : rule__ElkNode__Group_2__0__Impl rule__ElkNode__Group_2__1 ;
    public final void rule__ElkNode__Group_2__0() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:968:1: ( rule__ElkNode__Group_2__0__Impl rule__ElkNode__Group_2__1 )
            // InternalElkGraph.g:969:2: rule__ElkNode__Group_2__0__Impl rule__ElkNode__Group_2__1
            {
            pushFollow(FOLLOW_9);
            rule__ElkNode__Group_2__0__Impl();

            state._fsp--;

            pushFollow(FOLLOW_2);
            rule__ElkNode__Group_2__1();

            state._fsp--;


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ElkNode__Group_2__0"


    // $ANTLR start "rule__ElkNode__Group_2__0__Impl"
    // InternalElkGraph.g:976:1: rule__ElkNode__Group_2__0__Impl : ( '{' ) ;
    public final void rule__ElkNode__Group_2__0__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:980:1: ( ( '{' ) )
            // InternalElkGraph.g:981:1: ( '{' )
            {
            // InternalElkGraph.g:981:1: ( '{' )
            // InternalElkGraph.g:982:2: '{'
            {
             before(grammarAccess.getElkNodeAccess().getLeftCurlyBracketKeyword_2_0()); 
            match(input,18,FOLLOW_2); 
             after(grammarAccess.getElkNodeAccess().getLeftCurlyBracketKeyword_2_0()); 

            }


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ElkNode__Group_2__0__Impl"


    // $ANTLR start "rule__ElkNode__Group_2__1"
    // InternalElkGraph.g:991:1: rule__ElkNode__Group_2__1 : rule__ElkNode__Group_2__1__Impl rule__ElkNode__Group_2__2 ;
    public final void rule__ElkNode__Group_2__1() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:995:1: ( rule__ElkNode__Group_2__1__Impl rule__ElkNode__Group_2__2 )
            // InternalElkGraph.g:996:2: rule__ElkNode__Group_2__1__Impl rule__ElkNode__Group_2__2
            {
            pushFollow(FOLLOW_9);
            rule__ElkNode__Group_2__1__Impl();

            state._fsp--;

            pushFollow(FOLLOW_2);
            rule__ElkNode__Group_2__2();

            state._fsp--;


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ElkNode__Group_2__1"


    // $ANTLR start "rule__ElkNode__Group_2__1__Impl"
    // InternalElkGraph.g:1003:1: rule__ElkNode__Group_2__1__Impl : ( ( ruleShapeLayout )? ) ;
    public final void rule__ElkNode__Group_2__1__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:1007:1: ( ( ( ruleShapeLayout )? ) )
            // InternalElkGraph.g:1008:1: ( ( ruleShapeLayout )? )
            {
            // InternalElkGraph.g:1008:1: ( ( ruleShapeLayout )? )
            // InternalElkGraph.g:1009:2: ( ruleShapeLayout )?
            {
             before(grammarAccess.getElkNodeAccess().getShapeLayoutParserRuleCall_2_1()); 
            // InternalElkGraph.g:1010:2: ( ruleShapeLayout )?
            int alt14=2;
            int LA14_0 = input.LA(1);

            if ( (LA14_0==23) ) {
                alt14=1;
            }
            switch (alt14) {
                case 1 :
                    // InternalElkGraph.g:1010:3: ruleShapeLayout
                    {
                    pushFollow(FOLLOW_2);
                    ruleShapeLayout();

                    state._fsp--;


                    }
                    break;

            }

             after(grammarAccess.getElkNodeAccess().getShapeLayoutParserRuleCall_2_1()); 

            }


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ElkNode__Group_2__1__Impl"


    // $ANTLR start "rule__ElkNode__Group_2__2"
    // InternalElkGraph.g:1018:1: rule__ElkNode__Group_2__2 : rule__ElkNode__Group_2__2__Impl rule__ElkNode__Group_2__3 ;
    public final void rule__ElkNode__Group_2__2() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:1022:1: ( rule__ElkNode__Group_2__2__Impl rule__ElkNode__Group_2__3 )
            // InternalElkGraph.g:1023:2: rule__ElkNode__Group_2__2__Impl rule__ElkNode__Group_2__3
            {
            pushFollow(FOLLOW_9);
            rule__ElkNode__Group_2__2__Impl();

            state._fsp--;

            pushFollow(FOLLOW_2);
            rule__ElkNode__Group_2__3();

            state._fsp--;


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ElkNode__Group_2__2"


    // $ANTLR start "rule__ElkNode__Group_2__2__Impl"
    // InternalElkGraph.g:1030:1: rule__ElkNode__Group_2__2__Impl : ( ( rule__ElkNode__PropertiesAssignment_2_2 )* ) ;
    public final void rule__ElkNode__Group_2__2__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:1034:1: ( ( ( rule__ElkNode__PropertiesAssignment_2_2 )* ) )
            // InternalElkGraph.g:1035:1: ( ( rule__ElkNode__PropertiesAssignment_2_2 )* )
            {
            // InternalElkGraph.g:1035:1: ( ( rule__ElkNode__PropertiesAssignment_2_2 )* )
            // InternalElkGraph.g:1036:2: ( rule__ElkNode__PropertiesAssignment_2_2 )*
            {
             before(grammarAccess.getElkNodeAccess().getPropertiesAssignment_2_2()); 
            // InternalElkGraph.g:1037:2: ( rule__ElkNode__PropertiesAssignment_2_2 )*
            loop15:
            do {
                int alt15=2;
                int LA15_0 = input.LA(1);

                if ( (LA15_0==RULE_ID) ) {
                    alt15=1;
                }


                switch (alt15) {
            	case 1 :
            	    // InternalElkGraph.g:1037:3: rule__ElkNode__PropertiesAssignment_2_2
            	    {
            	    pushFollow(FOLLOW_5);
            	    rule__ElkNode__PropertiesAssignment_2_2();

            	    state._fsp--;


            	    }
            	    break;

            	default :
            	    break loop15;
                }
            } while (true);

             after(grammarAccess.getElkNodeAccess().getPropertiesAssignment_2_2()); 

            }


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ElkNode__Group_2__2__Impl"


    // $ANTLR start "rule__ElkNode__Group_2__3"
    // InternalElkGraph.g:1045:1: rule__ElkNode__Group_2__3 : rule__ElkNode__Group_2__3__Impl rule__ElkNode__Group_2__4 ;
    public final void rule__ElkNode__Group_2__3() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:1049:1: ( rule__ElkNode__Group_2__3__Impl rule__ElkNode__Group_2__4 )
            // InternalElkGraph.g:1050:2: rule__ElkNode__Group_2__3__Impl rule__ElkNode__Group_2__4
            {
            pushFollow(FOLLOW_9);
            rule__ElkNode__Group_2__3__Impl();

            state._fsp--;

            pushFollow(FOLLOW_2);
            rule__ElkNode__Group_2__4();

            state._fsp--;


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ElkNode__Group_2__3"


    // $ANTLR start "rule__ElkNode__Group_2__3__Impl"
    // InternalElkGraph.g:1057:1: rule__ElkNode__Group_2__3__Impl : ( ( rule__ElkNode__Alternatives_2_3 )* ) ;
    public final void rule__ElkNode__Group_2__3__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:1061:1: ( ( ( rule__ElkNode__Alternatives_2_3 )* ) )
            // InternalElkGraph.g:1062:1: ( ( rule__ElkNode__Alternatives_2_3 )* )
            {
            // InternalElkGraph.g:1062:1: ( ( rule__ElkNode__Alternatives_2_3 )* )
            // InternalElkGraph.g:1063:2: ( rule__ElkNode__Alternatives_2_3 )*
            {
             before(grammarAccess.getElkNodeAccess().getAlternatives_2_3()); 
            // InternalElkGraph.g:1064:2: ( rule__ElkNode__Alternatives_2_3 )*
            loop16:
            do {
                int alt16=2;
                int LA16_0 = input.LA(1);

                if ( (LA16_0==17||LA16_0==20||LA16_0==22||LA16_0==29) ) {
                    alt16=1;
                }


                switch (alt16) {
            	case 1 :
            	    // InternalElkGraph.g:1064:3: rule__ElkNode__Alternatives_2_3
            	    {
            	    pushFollow(FOLLOW_6);
            	    rule__ElkNode__Alternatives_2_3();

            	    state._fsp--;


            	    }
            	    break;

            	default :
            	    break loop16;
                }
            } while (true);

             after(grammarAccess.getElkNodeAccess().getAlternatives_2_3()); 

            }


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ElkNode__Group_2__3__Impl"


    // $ANTLR start "rule__ElkNode__Group_2__4"
    // InternalElkGraph.g:1072:1: rule__ElkNode__Group_2__4 : rule__ElkNode__Group_2__4__Impl ;
    public final void rule__ElkNode__Group_2__4() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:1076:1: ( rule__ElkNode__Group_2__4__Impl )
            // InternalElkGraph.g:1077:2: rule__ElkNode__Group_2__4__Impl
            {
            pushFollow(FOLLOW_2);
            rule__ElkNode__Group_2__4__Impl();

            state._fsp--;


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ElkNode__Group_2__4"


    // $ANTLR start "rule__ElkNode__Group_2__4__Impl"
    // InternalElkGraph.g:1083:1: rule__ElkNode__Group_2__4__Impl : ( '}' ) ;
    public final void rule__ElkNode__Group_2__4__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:1087:1: ( ( '}' ) )
            // InternalElkGraph.g:1088:1: ( '}' )
            {
            // InternalElkGraph.g:1088:1: ( '}' )
            // InternalElkGraph.g:1089:2: '}'
            {
             before(grammarAccess.getElkNodeAccess().getRightCurlyBracketKeyword_2_4()); 
            match(input,19,FOLLOW_2); 
             after(grammarAccess.getElkNodeAccess().getRightCurlyBracketKeyword_2_4()); 

            }


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ElkNode__Group_2__4__Impl"


    // $ANTLR start "rule__ElkLabel__Group__0"
    // InternalElkGraph.g:1099:1: rule__ElkLabel__Group__0 : rule__ElkLabel__Group__0__Impl rule__ElkLabel__Group__1 ;
    public final void rule__ElkLabel__Group__0() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:1103:1: ( rule__ElkLabel__Group__0__Impl rule__ElkLabel__Group__1 )
            // InternalElkGraph.g:1104:2: rule__ElkLabel__Group__0__Impl rule__ElkLabel__Group__1
            {
            pushFollow(FOLLOW_10);
            rule__ElkLabel__Group__0__Impl();

            state._fsp--;

            pushFollow(FOLLOW_2);
            rule__ElkLabel__Group__1();

            state._fsp--;


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ElkLabel__Group__0"


    // $ANTLR start "rule__ElkLabel__Group__0__Impl"
    // InternalElkGraph.g:1111:1: rule__ElkLabel__Group__0__Impl : ( 'label' ) ;
    public final void rule__ElkLabel__Group__0__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:1115:1: ( ( 'label' ) )
            // InternalElkGraph.g:1116:1: ( 'label' )
            {
            // InternalElkGraph.g:1116:1: ( 'label' )
            // InternalElkGraph.g:1117:2: 'label'
            {
             before(grammarAccess.getElkLabelAccess().getLabelKeyword_0()); 
            match(input,20,FOLLOW_2); 
             after(grammarAccess.getElkLabelAccess().getLabelKeyword_0()); 

            }


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ElkLabel__Group__0__Impl"


    // $ANTLR start "rule__ElkLabel__Group__1"
    // InternalElkGraph.g:1126:1: rule__ElkLabel__Group__1 : rule__ElkLabel__Group__1__Impl rule__ElkLabel__Group__2 ;
    public final void rule__ElkLabel__Group__1() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:1130:1: ( rule__ElkLabel__Group__1__Impl rule__ElkLabel__Group__2 )
            // InternalElkGraph.g:1131:2: rule__ElkLabel__Group__1__Impl rule__ElkLabel__Group__2
            {
            pushFollow(FOLLOW_10);
            rule__ElkLabel__Group__1__Impl();

            state._fsp--;

            pushFollow(FOLLOW_2);
            rule__ElkLabel__Group__2();

            state._fsp--;


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ElkLabel__Group__1"


    // $ANTLR start "rule__ElkLabel__Group__1__Impl"
    // InternalElkGraph.g:1138:1: rule__ElkLabel__Group__1__Impl : ( ( rule__ElkLabel__Group_1__0 )? ) ;
    public final void rule__ElkLabel__Group__1__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:1142:1: ( ( ( rule__ElkLabel__Group_1__0 )? ) )
            // InternalElkGraph.g:1143:1: ( ( rule__ElkLabel__Group_1__0 )? )
            {
            // InternalElkGraph.g:1143:1: ( ( rule__ElkLabel__Group_1__0 )? )
            // InternalElkGraph.g:1144:2: ( rule__ElkLabel__Group_1__0 )?
            {
             before(grammarAccess.getElkLabelAccess().getGroup_1()); 
            // InternalElkGraph.g:1145:2: ( rule__ElkLabel__Group_1__0 )?
            int alt17=2;
            int LA17_0 = input.LA(1);

            if ( (LA17_0==RULE_ID) ) {
                alt17=1;
            }
            switch (alt17) {
                case 1 :
                    // InternalElkGraph.g:1145:3: rule__ElkLabel__Group_1__0
                    {
                    pushFollow(FOLLOW_2);
                    rule__ElkLabel__Group_1__0();

                    state._fsp--;


                    }
                    break;

            }

             after(grammarAccess.getElkLabelAccess().getGroup_1()); 

            }


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ElkLabel__Group__1__Impl"


    // $ANTLR start "rule__ElkLabel__Group__2"
    // InternalElkGraph.g:1153:1: rule__ElkLabel__Group__2 : rule__ElkLabel__Group__2__Impl rule__ElkLabel__Group__3 ;
    public final void rule__ElkLabel__Group__2() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:1157:1: ( rule__ElkLabel__Group__2__Impl rule__ElkLabel__Group__3 )
            // InternalElkGraph.g:1158:2: rule__ElkLabel__Group__2__Impl rule__ElkLabel__Group__3
            {
            pushFollow(FOLLOW_8);
            rule__ElkLabel__Group__2__Impl();

            state._fsp--;

            pushFollow(FOLLOW_2);
            rule__ElkLabel__Group__3();

            state._fsp--;


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ElkLabel__Group__2"


    // $ANTLR start "rule__ElkLabel__Group__2__Impl"
    // InternalElkGraph.g:1165:1: rule__ElkLabel__Group__2__Impl : ( ( rule__ElkLabel__TextAssignment_2 ) ) ;
    public final void rule__ElkLabel__Group__2__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:1169:1: ( ( ( rule__ElkLabel__TextAssignment_2 ) ) )
            // InternalElkGraph.g:1170:1: ( ( rule__ElkLabel__TextAssignment_2 ) )
            {
            // InternalElkGraph.g:1170:1: ( ( rule__ElkLabel__TextAssignment_2 ) )
            // InternalElkGraph.g:1171:2: ( rule__ElkLabel__TextAssignment_2 )
            {
             before(grammarAccess.getElkLabelAccess().getTextAssignment_2()); 
            // InternalElkGraph.g:1172:2: ( rule__ElkLabel__TextAssignment_2 )
            // InternalElkGraph.g:1172:3: rule__ElkLabel__TextAssignment_2
            {
            pushFollow(FOLLOW_2);
            rule__ElkLabel__TextAssignment_2();

            state._fsp--;


            }

             after(grammarAccess.getElkLabelAccess().getTextAssignment_2()); 

            }


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ElkLabel__Group__2__Impl"


    // $ANTLR start "rule__ElkLabel__Group__3"
    // InternalElkGraph.g:1180:1: rule__ElkLabel__Group__3 : rule__ElkLabel__Group__3__Impl ;
    public final void rule__ElkLabel__Group__3() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:1184:1: ( rule__ElkLabel__Group__3__Impl )
            // InternalElkGraph.g:1185:2: rule__ElkLabel__Group__3__Impl
            {
            pushFollow(FOLLOW_2);
            rule__ElkLabel__Group__3__Impl();

            state._fsp--;


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ElkLabel__Group__3"


    // $ANTLR start "rule__ElkLabel__Group__3__Impl"
    // InternalElkGraph.g:1191:1: rule__ElkLabel__Group__3__Impl : ( ( rule__ElkLabel__Group_3__0 )? ) ;
    public final void rule__ElkLabel__Group__3__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:1195:1: ( ( ( rule__ElkLabel__Group_3__0 )? ) )
            // InternalElkGraph.g:1196:1: ( ( rule__ElkLabel__Group_3__0 )? )
            {
            // InternalElkGraph.g:1196:1: ( ( rule__ElkLabel__Group_3__0 )? )
            // InternalElkGraph.g:1197:2: ( rule__ElkLabel__Group_3__0 )?
            {
             before(grammarAccess.getElkLabelAccess().getGroup_3()); 
            // InternalElkGraph.g:1198:2: ( rule__ElkLabel__Group_3__0 )?
            int alt18=2;
            int LA18_0 = input.LA(1);

            if ( (LA18_0==18) ) {
                alt18=1;
            }
            switch (alt18) {
                case 1 :
                    // InternalElkGraph.g:1198:3: rule__ElkLabel__Group_3__0
                    {
                    pushFollow(FOLLOW_2);
                    rule__ElkLabel__Group_3__0();

                    state._fsp--;


                    }
                    break;

            }

             after(grammarAccess.getElkLabelAccess().getGroup_3()); 

            }


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ElkLabel__Group__3__Impl"


    // $ANTLR start "rule__ElkLabel__Group_1__0"
    // InternalElkGraph.g:1207:1: rule__ElkLabel__Group_1__0 : rule__ElkLabel__Group_1__0__Impl rule__ElkLabel__Group_1__1 ;
    public final void rule__ElkLabel__Group_1__0() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:1211:1: ( rule__ElkLabel__Group_1__0__Impl rule__ElkLabel__Group_1__1 )
            // InternalElkGraph.g:1212:2: rule__ElkLabel__Group_1__0__Impl rule__ElkLabel__Group_1__1
            {
            pushFollow(FOLLOW_11);
            rule__ElkLabel__Group_1__0__Impl();

            state._fsp--;

            pushFollow(FOLLOW_2);
            rule__ElkLabel__Group_1__1();

            state._fsp--;


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ElkLabel__Group_1__0"


    // $ANTLR start "rule__ElkLabel__Group_1__0__Impl"
    // InternalElkGraph.g:1219:1: rule__ElkLabel__Group_1__0__Impl : ( ( rule__ElkLabel__IdentifierAssignment_1_0 ) ) ;
    public final void rule__ElkLabel__Group_1__0__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:1223:1: ( ( ( rule__ElkLabel__IdentifierAssignment_1_0 ) ) )
            // InternalElkGraph.g:1224:1: ( ( rule__ElkLabel__IdentifierAssignment_1_0 ) )
            {
            // InternalElkGraph.g:1224:1: ( ( rule__ElkLabel__IdentifierAssignment_1_0 ) )
            // InternalElkGraph.g:1225:2: ( rule__ElkLabel__IdentifierAssignment_1_0 )
            {
             before(grammarAccess.getElkLabelAccess().getIdentifierAssignment_1_0()); 
            // InternalElkGraph.g:1226:2: ( rule__ElkLabel__IdentifierAssignment_1_0 )
            // InternalElkGraph.g:1226:3: rule__ElkLabel__IdentifierAssignment_1_0
            {
            pushFollow(FOLLOW_2);
            rule__ElkLabel__IdentifierAssignment_1_0();

            state._fsp--;


            }

             after(grammarAccess.getElkLabelAccess().getIdentifierAssignment_1_0()); 

            }


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ElkLabel__Group_1__0__Impl"


    // $ANTLR start "rule__ElkLabel__Group_1__1"
    // InternalElkGraph.g:1234:1: rule__ElkLabel__Group_1__1 : rule__ElkLabel__Group_1__1__Impl ;
    public final void rule__ElkLabel__Group_1__1() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:1238:1: ( rule__ElkLabel__Group_1__1__Impl )
            // InternalElkGraph.g:1239:2: rule__ElkLabel__Group_1__1__Impl
            {
            pushFollow(FOLLOW_2);
            rule__ElkLabel__Group_1__1__Impl();

            state._fsp--;


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ElkLabel__Group_1__1"


    // $ANTLR start "rule__ElkLabel__Group_1__1__Impl"
    // InternalElkGraph.g:1245:1: rule__ElkLabel__Group_1__1__Impl : ( ':' ) ;
    public final void rule__ElkLabel__Group_1__1__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:1249:1: ( ( ':' ) )
            // InternalElkGraph.g:1250:1: ( ':' )
            {
            // InternalElkGraph.g:1250:1: ( ':' )
            // InternalElkGraph.g:1251:2: ':'
            {
             before(grammarAccess.getElkLabelAccess().getColonKeyword_1_1()); 
            match(input,21,FOLLOW_2); 
             after(grammarAccess.getElkLabelAccess().getColonKeyword_1_1()); 

            }


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ElkLabel__Group_1__1__Impl"


    // $ANTLR start "rule__ElkLabel__Group_3__0"
    // InternalElkGraph.g:1261:1: rule__ElkLabel__Group_3__0 : rule__ElkLabel__Group_3__0__Impl rule__ElkLabel__Group_3__1 ;
    public final void rule__ElkLabel__Group_3__0() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:1265:1: ( rule__ElkLabel__Group_3__0__Impl rule__ElkLabel__Group_3__1 )
            // InternalElkGraph.g:1266:2: rule__ElkLabel__Group_3__0__Impl rule__ElkLabel__Group_3__1
            {
            pushFollow(FOLLOW_12);
            rule__ElkLabel__Group_3__0__Impl();

            state._fsp--;

            pushFollow(FOLLOW_2);
            rule__ElkLabel__Group_3__1();

            state._fsp--;


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ElkLabel__Group_3__0"


    // $ANTLR start "rule__ElkLabel__Group_3__0__Impl"
    // InternalElkGraph.g:1273:1: rule__ElkLabel__Group_3__0__Impl : ( '{' ) ;
    public final void rule__ElkLabel__Group_3__0__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:1277:1: ( ( '{' ) )
            // InternalElkGraph.g:1278:1: ( '{' )
            {
            // InternalElkGraph.g:1278:1: ( '{' )
            // InternalElkGraph.g:1279:2: '{'
            {
             before(grammarAccess.getElkLabelAccess().getLeftCurlyBracketKeyword_3_0()); 
            match(input,18,FOLLOW_2); 
             after(grammarAccess.getElkLabelAccess().getLeftCurlyBracketKeyword_3_0()); 

            }


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ElkLabel__Group_3__0__Impl"


    // $ANTLR start "rule__ElkLabel__Group_3__1"
    // InternalElkGraph.g:1288:1: rule__ElkLabel__Group_3__1 : rule__ElkLabel__Group_3__1__Impl rule__ElkLabel__Group_3__2 ;
    public final void rule__ElkLabel__Group_3__1() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:1292:1: ( rule__ElkLabel__Group_3__1__Impl rule__ElkLabel__Group_3__2 )
            // InternalElkGraph.g:1293:2: rule__ElkLabel__Group_3__1__Impl rule__ElkLabel__Group_3__2
            {
            pushFollow(FOLLOW_12);
            rule__ElkLabel__Group_3__1__Impl();

            state._fsp--;

            pushFollow(FOLLOW_2);
            rule__ElkLabel__Group_3__2();

            state._fsp--;


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ElkLabel__Group_3__1"


    // $ANTLR start "rule__ElkLabel__Group_3__1__Impl"
    // InternalElkGraph.g:1300:1: rule__ElkLabel__Group_3__1__Impl : ( ( ruleShapeLayout )? ) ;
    public final void rule__ElkLabel__Group_3__1__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:1304:1: ( ( ( ruleShapeLayout )? ) )
            // InternalElkGraph.g:1305:1: ( ( ruleShapeLayout )? )
            {
            // InternalElkGraph.g:1305:1: ( ( ruleShapeLayout )? )
            // InternalElkGraph.g:1306:2: ( ruleShapeLayout )?
            {
             before(grammarAccess.getElkLabelAccess().getShapeLayoutParserRuleCall_3_1()); 
            // InternalElkGraph.g:1307:2: ( ruleShapeLayout )?
            int alt19=2;
            int LA19_0 = input.LA(1);

            if ( (LA19_0==23) ) {
                alt19=1;
            }
            switch (alt19) {
                case 1 :
                    // InternalElkGraph.g:1307:3: ruleShapeLayout
                    {
                    pushFollow(FOLLOW_2);
                    ruleShapeLayout();

                    state._fsp--;


                    }
                    break;

            }

             after(grammarAccess.getElkLabelAccess().getShapeLayoutParserRuleCall_3_1()); 

            }


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ElkLabel__Group_3__1__Impl"


    // $ANTLR start "rule__ElkLabel__Group_3__2"
    // InternalElkGraph.g:1315:1: rule__ElkLabel__Group_3__2 : rule__ElkLabel__Group_3__2__Impl rule__ElkLabel__Group_3__3 ;
    public final void rule__ElkLabel__Group_3__2() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:1319:1: ( rule__ElkLabel__Group_3__2__Impl rule__ElkLabel__Group_3__3 )
            // InternalElkGraph.g:1320:2: rule__ElkLabel__Group_3__2__Impl rule__ElkLabel__Group_3__3
            {
            pushFollow(FOLLOW_12);
            rule__ElkLabel__Group_3__2__Impl();

            state._fsp--;

            pushFollow(FOLLOW_2);
            rule__ElkLabel__Group_3__3();

            state._fsp--;


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ElkLabel__Group_3__2"


    // $ANTLR start "rule__ElkLabel__Group_3__2__Impl"
    // InternalElkGraph.g:1327:1: rule__ElkLabel__Group_3__2__Impl : ( ( rule__ElkLabel__PropertiesAssignment_3_2 )* ) ;
    public final void rule__ElkLabel__Group_3__2__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:1331:1: ( ( ( rule__ElkLabel__PropertiesAssignment_3_2 )* ) )
            // InternalElkGraph.g:1332:1: ( ( rule__ElkLabel__PropertiesAssignment_3_2 )* )
            {
            // InternalElkGraph.g:1332:1: ( ( rule__ElkLabel__PropertiesAssignment_3_2 )* )
            // InternalElkGraph.g:1333:2: ( rule__ElkLabel__PropertiesAssignment_3_2 )*
            {
             before(grammarAccess.getElkLabelAccess().getPropertiesAssignment_3_2()); 
            // InternalElkGraph.g:1334:2: ( rule__ElkLabel__PropertiesAssignment_3_2 )*
            loop20:
            do {
                int alt20=2;
                int LA20_0 = input.LA(1);

                if ( (LA20_0==RULE_ID) ) {
                    alt20=1;
                }


                switch (alt20) {
            	case 1 :
            	    // InternalElkGraph.g:1334:3: rule__ElkLabel__PropertiesAssignment_3_2
            	    {
            	    pushFollow(FOLLOW_5);
            	    rule__ElkLabel__PropertiesAssignment_3_2();

            	    state._fsp--;


            	    }
            	    break;

            	default :
            	    break loop20;
                }
            } while (true);

             after(grammarAccess.getElkLabelAccess().getPropertiesAssignment_3_2()); 

            }


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ElkLabel__Group_3__2__Impl"


    // $ANTLR start "rule__ElkLabel__Group_3__3"
    // InternalElkGraph.g:1342:1: rule__ElkLabel__Group_3__3 : rule__ElkLabel__Group_3__3__Impl rule__ElkLabel__Group_3__4 ;
    public final void rule__ElkLabel__Group_3__3() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:1346:1: ( rule__ElkLabel__Group_3__3__Impl rule__ElkLabel__Group_3__4 )
            // InternalElkGraph.g:1347:2: rule__ElkLabel__Group_3__3__Impl rule__ElkLabel__Group_3__4
            {
            pushFollow(FOLLOW_12);
            rule__ElkLabel__Group_3__3__Impl();

            state._fsp--;

            pushFollow(FOLLOW_2);
            rule__ElkLabel__Group_3__4();

            state._fsp--;


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ElkLabel__Group_3__3"


    // $ANTLR start "rule__ElkLabel__Group_3__3__Impl"
    // InternalElkGraph.g:1354:1: rule__ElkLabel__Group_3__3__Impl : ( ( rule__ElkLabel__LabelsAssignment_3_3 )* ) ;
    public final void rule__ElkLabel__Group_3__3__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:1358:1: ( ( ( rule__ElkLabel__LabelsAssignment_3_3 )* ) )
            // InternalElkGraph.g:1359:1: ( ( rule__ElkLabel__LabelsAssignment_3_3 )* )
            {
            // InternalElkGraph.g:1359:1: ( ( rule__ElkLabel__LabelsAssignment_3_3 )* )
            // InternalElkGraph.g:1360:2: ( rule__ElkLabel__LabelsAssignment_3_3 )*
            {
             before(grammarAccess.getElkLabelAccess().getLabelsAssignment_3_3()); 
            // InternalElkGraph.g:1361:2: ( rule__ElkLabel__LabelsAssignment_3_3 )*
            loop21:
            do {
                int alt21=2;
                int LA21_0 = input.LA(1);

                if ( (LA21_0==20) ) {
                    alt21=1;
                }


                switch (alt21) {
            	case 1 :
            	    // InternalElkGraph.g:1361:3: rule__ElkLabel__LabelsAssignment_3_3
            	    {
            	    pushFollow(FOLLOW_13);
            	    rule__ElkLabel__LabelsAssignment_3_3();

            	    state._fsp--;


            	    }
            	    break;

            	default :
            	    break loop21;
                }
            } while (true);

             after(grammarAccess.getElkLabelAccess().getLabelsAssignment_3_3()); 

            }


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ElkLabel__Group_3__3__Impl"


    // $ANTLR start "rule__ElkLabel__Group_3__4"
    // InternalElkGraph.g:1369:1: rule__ElkLabel__Group_3__4 : rule__ElkLabel__Group_3__4__Impl ;
    public final void rule__ElkLabel__Group_3__4() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:1373:1: ( rule__ElkLabel__Group_3__4__Impl )
            // InternalElkGraph.g:1374:2: rule__ElkLabel__Group_3__4__Impl
            {
            pushFollow(FOLLOW_2);
            rule__ElkLabel__Group_3__4__Impl();

            state._fsp--;


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ElkLabel__Group_3__4"


    // $ANTLR start "rule__ElkLabel__Group_3__4__Impl"
    // InternalElkGraph.g:1380:1: rule__ElkLabel__Group_3__4__Impl : ( '}' ) ;
    public final void rule__ElkLabel__Group_3__4__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:1384:1: ( ( '}' ) )
            // InternalElkGraph.g:1385:1: ( '}' )
            {
            // InternalElkGraph.g:1385:1: ( '}' )
            // InternalElkGraph.g:1386:2: '}'
            {
             before(grammarAccess.getElkLabelAccess().getRightCurlyBracketKeyword_3_4()); 
            match(input,19,FOLLOW_2); 
             after(grammarAccess.getElkLabelAccess().getRightCurlyBracketKeyword_3_4()); 

            }


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ElkLabel__Group_3__4__Impl"


    // $ANTLR start "rule__ElkPort__Group__0"
    // InternalElkGraph.g:1396:1: rule__ElkPort__Group__0 : rule__ElkPort__Group__0__Impl rule__ElkPort__Group__1 ;
    public final void rule__ElkPort__Group__0() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:1400:1: ( rule__ElkPort__Group__0__Impl rule__ElkPort__Group__1 )
            // InternalElkGraph.g:1401:2: rule__ElkPort__Group__0__Impl rule__ElkPort__Group__1
            {
            pushFollow(FOLLOW_7);
            rule__ElkPort__Group__0__Impl();

            state._fsp--;

            pushFollow(FOLLOW_2);
            rule__ElkPort__Group__1();

            state._fsp--;


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ElkPort__Group__0"


    // $ANTLR start "rule__ElkPort__Group__0__Impl"
    // InternalElkGraph.g:1408:1: rule__ElkPort__Group__0__Impl : ( 'port' ) ;
    public final void rule__ElkPort__Group__0__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:1412:1: ( ( 'port' ) )
            // InternalElkGraph.g:1413:1: ( 'port' )
            {
            // InternalElkGraph.g:1413:1: ( 'port' )
            // InternalElkGraph.g:1414:2: 'port'
            {
             before(grammarAccess.getElkPortAccess().getPortKeyword_0()); 
            match(input,22,FOLLOW_2); 
             after(grammarAccess.getElkPortAccess().getPortKeyword_0()); 

            }


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ElkPort__Group__0__Impl"


    // $ANTLR start "rule__ElkPort__Group__1"
    // InternalElkGraph.g:1423:1: rule__ElkPort__Group__1 : rule__ElkPort__Group__1__Impl rule__ElkPort__Group__2 ;
    public final void rule__ElkPort__Group__1() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:1427:1: ( rule__ElkPort__Group__1__Impl rule__ElkPort__Group__2 )
            // InternalElkGraph.g:1428:2: rule__ElkPort__Group__1__Impl rule__ElkPort__Group__2
            {
            pushFollow(FOLLOW_8);
            rule__ElkPort__Group__1__Impl();

            state._fsp--;

            pushFollow(FOLLOW_2);
            rule__ElkPort__Group__2();

            state._fsp--;


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ElkPort__Group__1"


    // $ANTLR start "rule__ElkPort__Group__1__Impl"
    // InternalElkGraph.g:1435:1: rule__ElkPort__Group__1__Impl : ( ( rule__ElkPort__IdentifierAssignment_1 ) ) ;
    public final void rule__ElkPort__Group__1__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:1439:1: ( ( ( rule__ElkPort__IdentifierAssignment_1 ) ) )
            // InternalElkGraph.g:1440:1: ( ( rule__ElkPort__IdentifierAssignment_1 ) )
            {
            // InternalElkGraph.g:1440:1: ( ( rule__ElkPort__IdentifierAssignment_1 ) )
            // InternalElkGraph.g:1441:2: ( rule__ElkPort__IdentifierAssignment_1 )
            {
             before(grammarAccess.getElkPortAccess().getIdentifierAssignment_1()); 
            // InternalElkGraph.g:1442:2: ( rule__ElkPort__IdentifierAssignment_1 )
            // InternalElkGraph.g:1442:3: rule__ElkPort__IdentifierAssignment_1
            {
            pushFollow(FOLLOW_2);
            rule__ElkPort__IdentifierAssignment_1();

            state._fsp--;


            }

             after(grammarAccess.getElkPortAccess().getIdentifierAssignment_1()); 

            }


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ElkPort__Group__1__Impl"


    // $ANTLR start "rule__ElkPort__Group__2"
    // InternalElkGraph.g:1450:1: rule__ElkPort__Group__2 : rule__ElkPort__Group__2__Impl ;
    public final void rule__ElkPort__Group__2() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:1454:1: ( rule__ElkPort__Group__2__Impl )
            // InternalElkGraph.g:1455:2: rule__ElkPort__Group__2__Impl
            {
            pushFollow(FOLLOW_2);
            rule__ElkPort__Group__2__Impl();

            state._fsp--;


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ElkPort__Group__2"


    // $ANTLR start "rule__ElkPort__Group__2__Impl"
    // InternalElkGraph.g:1461:1: rule__ElkPort__Group__2__Impl : ( ( rule__ElkPort__Group_2__0 )? ) ;
    public final void rule__ElkPort__Group__2__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:1465:1: ( ( ( rule__ElkPort__Group_2__0 )? ) )
            // InternalElkGraph.g:1466:1: ( ( rule__ElkPort__Group_2__0 )? )
            {
            // InternalElkGraph.g:1466:1: ( ( rule__ElkPort__Group_2__0 )? )
            // InternalElkGraph.g:1467:2: ( rule__ElkPort__Group_2__0 )?
            {
             before(grammarAccess.getElkPortAccess().getGroup_2()); 
            // InternalElkGraph.g:1468:2: ( rule__ElkPort__Group_2__0 )?
            int alt22=2;
            int LA22_0 = input.LA(1);

            if ( (LA22_0==18) ) {
                alt22=1;
            }
            switch (alt22) {
                case 1 :
                    // InternalElkGraph.g:1468:3: rule__ElkPort__Group_2__0
                    {
                    pushFollow(FOLLOW_2);
                    rule__ElkPort__Group_2__0();

                    state._fsp--;


                    }
                    break;

            }

             after(grammarAccess.getElkPortAccess().getGroup_2()); 

            }


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ElkPort__Group__2__Impl"


    // $ANTLR start "rule__ElkPort__Group_2__0"
    // InternalElkGraph.g:1477:1: rule__ElkPort__Group_2__0 : rule__ElkPort__Group_2__0__Impl rule__ElkPort__Group_2__1 ;
    public final void rule__ElkPort__Group_2__0() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:1481:1: ( rule__ElkPort__Group_2__0__Impl rule__ElkPort__Group_2__1 )
            // InternalElkGraph.g:1482:2: rule__ElkPort__Group_2__0__Impl rule__ElkPort__Group_2__1
            {
            pushFollow(FOLLOW_12);
            rule__ElkPort__Group_2__0__Impl();

            state._fsp--;

            pushFollow(FOLLOW_2);
            rule__ElkPort__Group_2__1();

            state._fsp--;


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ElkPort__Group_2__0"


    // $ANTLR start "rule__ElkPort__Group_2__0__Impl"
    // InternalElkGraph.g:1489:1: rule__ElkPort__Group_2__0__Impl : ( '{' ) ;
    public final void rule__ElkPort__Group_2__0__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:1493:1: ( ( '{' ) )
            // InternalElkGraph.g:1494:1: ( '{' )
            {
            // InternalElkGraph.g:1494:1: ( '{' )
            // InternalElkGraph.g:1495:2: '{'
            {
             before(grammarAccess.getElkPortAccess().getLeftCurlyBracketKeyword_2_0()); 
            match(input,18,FOLLOW_2); 
             after(grammarAccess.getElkPortAccess().getLeftCurlyBracketKeyword_2_0()); 

            }


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ElkPort__Group_2__0__Impl"


    // $ANTLR start "rule__ElkPort__Group_2__1"
    // InternalElkGraph.g:1504:1: rule__ElkPort__Group_2__1 : rule__ElkPort__Group_2__1__Impl rule__ElkPort__Group_2__2 ;
    public final void rule__ElkPort__Group_2__1() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:1508:1: ( rule__ElkPort__Group_2__1__Impl rule__ElkPort__Group_2__2 )
            // InternalElkGraph.g:1509:2: rule__ElkPort__Group_2__1__Impl rule__ElkPort__Group_2__2
            {
            pushFollow(FOLLOW_12);
            rule__ElkPort__Group_2__1__Impl();

            state._fsp--;

            pushFollow(FOLLOW_2);
            rule__ElkPort__Group_2__2();

            state._fsp--;


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ElkPort__Group_2__1"


    // $ANTLR start "rule__ElkPort__Group_2__1__Impl"
    // InternalElkGraph.g:1516:1: rule__ElkPort__Group_2__1__Impl : ( ( ruleShapeLayout )? ) ;
    public final void rule__ElkPort__Group_2__1__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:1520:1: ( ( ( ruleShapeLayout )? ) )
            // InternalElkGraph.g:1521:1: ( ( ruleShapeLayout )? )
            {
            // InternalElkGraph.g:1521:1: ( ( ruleShapeLayout )? )
            // InternalElkGraph.g:1522:2: ( ruleShapeLayout )?
            {
             before(grammarAccess.getElkPortAccess().getShapeLayoutParserRuleCall_2_1()); 
            // InternalElkGraph.g:1523:2: ( ruleShapeLayout )?
            int alt23=2;
            int LA23_0 = input.LA(1);

            if ( (LA23_0==23) ) {
                alt23=1;
            }
            switch (alt23) {
                case 1 :
                    // InternalElkGraph.g:1523:3: ruleShapeLayout
                    {
                    pushFollow(FOLLOW_2);
                    ruleShapeLayout();

                    state._fsp--;


                    }
                    break;

            }

             after(grammarAccess.getElkPortAccess().getShapeLayoutParserRuleCall_2_1()); 

            }


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ElkPort__Group_2__1__Impl"


    // $ANTLR start "rule__ElkPort__Group_2__2"
    // InternalElkGraph.g:1531:1: rule__ElkPort__Group_2__2 : rule__ElkPort__Group_2__2__Impl rule__ElkPort__Group_2__3 ;
    public final void rule__ElkPort__Group_2__2() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:1535:1: ( rule__ElkPort__Group_2__2__Impl rule__ElkPort__Group_2__3 )
            // InternalElkGraph.g:1536:2: rule__ElkPort__Group_2__2__Impl rule__ElkPort__Group_2__3
            {
            pushFollow(FOLLOW_12);
            rule__ElkPort__Group_2__2__Impl();

            state._fsp--;

            pushFollow(FOLLOW_2);
            rule__ElkPort__Group_2__3();

            state._fsp--;


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ElkPort__Group_2__2"


    // $ANTLR start "rule__ElkPort__Group_2__2__Impl"
    // InternalElkGraph.g:1543:1: rule__ElkPort__Group_2__2__Impl : ( ( rule__ElkPort__PropertiesAssignment_2_2 )* ) ;
    public final void rule__ElkPort__Group_2__2__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:1547:1: ( ( ( rule__ElkPort__PropertiesAssignment_2_2 )* ) )
            // InternalElkGraph.g:1548:1: ( ( rule__ElkPort__PropertiesAssignment_2_2 )* )
            {
            // InternalElkGraph.g:1548:1: ( ( rule__ElkPort__PropertiesAssignment_2_2 )* )
            // InternalElkGraph.g:1549:2: ( rule__ElkPort__PropertiesAssignment_2_2 )*
            {
             before(grammarAccess.getElkPortAccess().getPropertiesAssignment_2_2()); 
            // InternalElkGraph.g:1550:2: ( rule__ElkPort__PropertiesAssignment_2_2 )*
            loop24:
            do {
                int alt24=2;
                int LA24_0 = input.LA(1);

                if ( (LA24_0==RULE_ID) ) {
                    alt24=1;
                }


                switch (alt24) {
            	case 1 :
            	    // InternalElkGraph.g:1550:3: rule__ElkPort__PropertiesAssignment_2_2
            	    {
            	    pushFollow(FOLLOW_5);
            	    rule__ElkPort__PropertiesAssignment_2_2();

            	    state._fsp--;


            	    }
            	    break;

            	default :
            	    break loop24;
                }
            } while (true);

             after(grammarAccess.getElkPortAccess().getPropertiesAssignment_2_2()); 

            }


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ElkPort__Group_2__2__Impl"


    // $ANTLR start "rule__ElkPort__Group_2__3"
    // InternalElkGraph.g:1558:1: rule__ElkPort__Group_2__3 : rule__ElkPort__Group_2__3__Impl rule__ElkPort__Group_2__4 ;
    public final void rule__ElkPort__Group_2__3() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:1562:1: ( rule__ElkPort__Group_2__3__Impl rule__ElkPort__Group_2__4 )
            // InternalElkGraph.g:1563:2: rule__ElkPort__Group_2__3__Impl rule__ElkPort__Group_2__4
            {
            pushFollow(FOLLOW_12);
            rule__ElkPort__Group_2__3__Impl();

            state._fsp--;

            pushFollow(FOLLOW_2);
            rule__ElkPort__Group_2__4();

            state._fsp--;


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ElkPort__Group_2__3"


    // $ANTLR start "rule__ElkPort__Group_2__3__Impl"
    // InternalElkGraph.g:1570:1: rule__ElkPort__Group_2__3__Impl : ( ( rule__ElkPort__LabelsAssignment_2_3 )* ) ;
    public final void rule__ElkPort__Group_2__3__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:1574:1: ( ( ( rule__ElkPort__LabelsAssignment_2_3 )* ) )
            // InternalElkGraph.g:1575:1: ( ( rule__ElkPort__LabelsAssignment_2_3 )* )
            {
            // InternalElkGraph.g:1575:1: ( ( rule__ElkPort__LabelsAssignment_2_3 )* )
            // InternalElkGraph.g:1576:2: ( rule__ElkPort__LabelsAssignment_2_3 )*
            {
             before(grammarAccess.getElkPortAccess().getLabelsAssignment_2_3()); 
            // InternalElkGraph.g:1577:2: ( rule__ElkPort__LabelsAssignment_2_3 )*
            loop25:
            do {
                int alt25=2;
                int LA25_0 = input.LA(1);

                if ( (LA25_0==20) ) {
                    alt25=1;
                }


                switch (alt25) {
            	case 1 :
            	    // InternalElkGraph.g:1577:3: rule__ElkPort__LabelsAssignment_2_3
            	    {
            	    pushFollow(FOLLOW_13);
            	    rule__ElkPort__LabelsAssignment_2_3();

            	    state._fsp--;


            	    }
            	    break;

            	default :
            	    break loop25;
                }
            } while (true);

             after(grammarAccess.getElkPortAccess().getLabelsAssignment_2_3()); 

            }


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ElkPort__Group_2__3__Impl"


    // $ANTLR start "rule__ElkPort__Group_2__4"
    // InternalElkGraph.g:1585:1: rule__ElkPort__Group_2__4 : rule__ElkPort__Group_2__4__Impl ;
    public final void rule__ElkPort__Group_2__4() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:1589:1: ( rule__ElkPort__Group_2__4__Impl )
            // InternalElkGraph.g:1590:2: rule__ElkPort__Group_2__4__Impl
            {
            pushFollow(FOLLOW_2);
            rule__ElkPort__Group_2__4__Impl();

            state._fsp--;


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ElkPort__Group_2__4"


    // $ANTLR start "rule__ElkPort__Group_2__4__Impl"
    // InternalElkGraph.g:1596:1: rule__ElkPort__Group_2__4__Impl : ( '}' ) ;
    public final void rule__ElkPort__Group_2__4__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:1600:1: ( ( '}' ) )
            // InternalElkGraph.g:1601:1: ( '}' )
            {
            // InternalElkGraph.g:1601:1: ( '}' )
            // InternalElkGraph.g:1602:2: '}'
            {
             before(grammarAccess.getElkPortAccess().getRightCurlyBracketKeyword_2_4()); 
            match(input,19,FOLLOW_2); 
             after(grammarAccess.getElkPortAccess().getRightCurlyBracketKeyword_2_4()); 

            }


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ElkPort__Group_2__4__Impl"


    // $ANTLR start "rule__ShapeLayout__Group__0"
    // InternalElkGraph.g:1612:1: rule__ShapeLayout__Group__0 : rule__ShapeLayout__Group__0__Impl rule__ShapeLayout__Group__1 ;
    public final void rule__ShapeLayout__Group__0() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:1616:1: ( rule__ShapeLayout__Group__0__Impl rule__ShapeLayout__Group__1 )
            // InternalElkGraph.g:1617:2: rule__ShapeLayout__Group__0__Impl rule__ShapeLayout__Group__1
            {
            pushFollow(FOLLOW_14);
            rule__ShapeLayout__Group__0__Impl();

            state._fsp--;

            pushFollow(FOLLOW_2);
            rule__ShapeLayout__Group__1();

            state._fsp--;


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ShapeLayout__Group__0"


    // $ANTLR start "rule__ShapeLayout__Group__0__Impl"
    // InternalElkGraph.g:1624:1: rule__ShapeLayout__Group__0__Impl : ( 'layout' ) ;
    public final void rule__ShapeLayout__Group__0__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:1628:1: ( ( 'layout' ) )
            // InternalElkGraph.g:1629:1: ( 'layout' )
            {
            // InternalElkGraph.g:1629:1: ( 'layout' )
            // InternalElkGraph.g:1630:2: 'layout'
            {
             before(grammarAccess.getShapeLayoutAccess().getLayoutKeyword_0()); 
            match(input,23,FOLLOW_2); 
             after(grammarAccess.getShapeLayoutAccess().getLayoutKeyword_0()); 

            }


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ShapeLayout__Group__0__Impl"


    // $ANTLR start "rule__ShapeLayout__Group__1"
    // InternalElkGraph.g:1639:1: rule__ShapeLayout__Group__1 : rule__ShapeLayout__Group__1__Impl rule__ShapeLayout__Group__2 ;
    public final void rule__ShapeLayout__Group__1() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:1643:1: ( rule__ShapeLayout__Group__1__Impl rule__ShapeLayout__Group__2 )
            // InternalElkGraph.g:1644:2: rule__ShapeLayout__Group__1__Impl rule__ShapeLayout__Group__2
            {
            pushFollow(FOLLOW_15);
            rule__ShapeLayout__Group__1__Impl();

            state._fsp--;

            pushFollow(FOLLOW_2);
            rule__ShapeLayout__Group__2();

            state._fsp--;


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ShapeLayout__Group__1"


    // $ANTLR start "rule__ShapeLayout__Group__1__Impl"
    // InternalElkGraph.g:1651:1: rule__ShapeLayout__Group__1__Impl : ( '[' ) ;
    public final void rule__ShapeLayout__Group__1__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:1655:1: ( ( '[' ) )
            // InternalElkGraph.g:1656:1: ( '[' )
            {
            // InternalElkGraph.g:1656:1: ( '[' )
            // InternalElkGraph.g:1657:2: '['
            {
             before(grammarAccess.getShapeLayoutAccess().getLeftSquareBracketKeyword_1()); 
            match(input,24,FOLLOW_2); 
             after(grammarAccess.getShapeLayoutAccess().getLeftSquareBracketKeyword_1()); 

            }


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ShapeLayout__Group__1__Impl"


    // $ANTLR start "rule__ShapeLayout__Group__2"
    // InternalElkGraph.g:1666:1: rule__ShapeLayout__Group__2 : rule__ShapeLayout__Group__2__Impl rule__ShapeLayout__Group__3 ;
    public final void rule__ShapeLayout__Group__2() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:1670:1: ( rule__ShapeLayout__Group__2__Impl rule__ShapeLayout__Group__3 )
            // InternalElkGraph.g:1671:2: rule__ShapeLayout__Group__2__Impl rule__ShapeLayout__Group__3
            {
            pushFollow(FOLLOW_16);
            rule__ShapeLayout__Group__2__Impl();

            state._fsp--;

            pushFollow(FOLLOW_2);
            rule__ShapeLayout__Group__3();

            state._fsp--;


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ShapeLayout__Group__2"


    // $ANTLR start "rule__ShapeLayout__Group__2__Impl"
    // InternalElkGraph.g:1678:1: rule__ShapeLayout__Group__2__Impl : ( ( rule__ShapeLayout__UnorderedGroup_2 ) ) ;
    public final void rule__ShapeLayout__Group__2__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:1682:1: ( ( ( rule__ShapeLayout__UnorderedGroup_2 ) ) )
            // InternalElkGraph.g:1683:1: ( ( rule__ShapeLayout__UnorderedGroup_2 ) )
            {
            // InternalElkGraph.g:1683:1: ( ( rule__ShapeLayout__UnorderedGroup_2 ) )
            // InternalElkGraph.g:1684:2: ( rule__ShapeLayout__UnorderedGroup_2 )
            {
             before(grammarAccess.getShapeLayoutAccess().getUnorderedGroup_2()); 
            // InternalElkGraph.g:1685:2: ( rule__ShapeLayout__UnorderedGroup_2 )
            // InternalElkGraph.g:1685:3: rule__ShapeLayout__UnorderedGroup_2
            {
            pushFollow(FOLLOW_2);
            rule__ShapeLayout__UnorderedGroup_2();

            state._fsp--;


            }

             after(grammarAccess.getShapeLayoutAccess().getUnorderedGroup_2()); 

            }


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ShapeLayout__Group__2__Impl"


    // $ANTLR start "rule__ShapeLayout__Group__3"
    // InternalElkGraph.g:1693:1: rule__ShapeLayout__Group__3 : rule__ShapeLayout__Group__3__Impl ;
    public final void rule__ShapeLayout__Group__3() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:1697:1: ( rule__ShapeLayout__Group__3__Impl )
            // InternalElkGraph.g:1698:2: rule__ShapeLayout__Group__3__Impl
            {
            pushFollow(FOLLOW_2);
            rule__ShapeLayout__Group__3__Impl();

            state._fsp--;


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ShapeLayout__Group__3"


    // $ANTLR start "rule__ShapeLayout__Group__3__Impl"
    // InternalElkGraph.g:1704:1: rule__ShapeLayout__Group__3__Impl : ( ']' ) ;
    public final void rule__ShapeLayout__Group__3__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:1708:1: ( ( ']' ) )
            // InternalElkGraph.g:1709:1: ( ']' )
            {
            // InternalElkGraph.g:1709:1: ( ']' )
            // InternalElkGraph.g:1710:2: ']'
            {
             before(grammarAccess.getShapeLayoutAccess().getRightSquareBracketKeyword_3()); 
            match(input,25,FOLLOW_2); 
             after(grammarAccess.getShapeLayoutAccess().getRightSquareBracketKeyword_3()); 

            }


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ShapeLayout__Group__3__Impl"


    // $ANTLR start "rule__ShapeLayout__Group_2_0__0"
    // InternalElkGraph.g:1720:1: rule__ShapeLayout__Group_2_0__0 : rule__ShapeLayout__Group_2_0__0__Impl rule__ShapeLayout__Group_2_0__1 ;
    public final void rule__ShapeLayout__Group_2_0__0() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:1724:1: ( rule__ShapeLayout__Group_2_0__0__Impl rule__ShapeLayout__Group_2_0__1 )
            // InternalElkGraph.g:1725:2: rule__ShapeLayout__Group_2_0__0__Impl rule__ShapeLayout__Group_2_0__1
            {
            pushFollow(FOLLOW_11);
            rule__ShapeLayout__Group_2_0__0__Impl();

            state._fsp--;

            pushFollow(FOLLOW_2);
            rule__ShapeLayout__Group_2_0__1();

            state._fsp--;


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ShapeLayout__Group_2_0__0"


    // $ANTLR start "rule__ShapeLayout__Group_2_0__0__Impl"
    // InternalElkGraph.g:1732:1: rule__ShapeLayout__Group_2_0__0__Impl : ( 'position' ) ;
    public final void rule__ShapeLayout__Group_2_0__0__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:1736:1: ( ( 'position' ) )
            // InternalElkGraph.g:1737:1: ( 'position' )
            {
            // InternalElkGraph.g:1737:1: ( 'position' )
            // InternalElkGraph.g:1738:2: 'position'
            {
             before(grammarAccess.getShapeLayoutAccess().getPositionKeyword_2_0_0()); 
            match(input,26,FOLLOW_2); 
             after(grammarAccess.getShapeLayoutAccess().getPositionKeyword_2_0_0()); 

            }


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ShapeLayout__Group_2_0__0__Impl"


    // $ANTLR start "rule__ShapeLayout__Group_2_0__1"
    // InternalElkGraph.g:1747:1: rule__ShapeLayout__Group_2_0__1 : rule__ShapeLayout__Group_2_0__1__Impl rule__ShapeLayout__Group_2_0__2 ;
    public final void rule__ShapeLayout__Group_2_0__1() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:1751:1: ( rule__ShapeLayout__Group_2_0__1__Impl rule__ShapeLayout__Group_2_0__2 )
            // InternalElkGraph.g:1752:2: rule__ShapeLayout__Group_2_0__1__Impl rule__ShapeLayout__Group_2_0__2
            {
            pushFollow(FOLLOW_17);
            rule__ShapeLayout__Group_2_0__1__Impl();

            state._fsp--;

            pushFollow(FOLLOW_2);
            rule__ShapeLayout__Group_2_0__2();

            state._fsp--;


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ShapeLayout__Group_2_0__1"


    // $ANTLR start "rule__ShapeLayout__Group_2_0__1__Impl"
    // InternalElkGraph.g:1759:1: rule__ShapeLayout__Group_2_0__1__Impl : ( ':' ) ;
    public final void rule__ShapeLayout__Group_2_0__1__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:1763:1: ( ( ':' ) )
            // InternalElkGraph.g:1764:1: ( ':' )
            {
            // InternalElkGraph.g:1764:1: ( ':' )
            // InternalElkGraph.g:1765:2: ':'
            {
             before(grammarAccess.getShapeLayoutAccess().getColonKeyword_2_0_1()); 
            match(input,21,FOLLOW_2); 
             after(grammarAccess.getShapeLayoutAccess().getColonKeyword_2_0_1()); 

            }


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ShapeLayout__Group_2_0__1__Impl"


    // $ANTLR start "rule__ShapeLayout__Group_2_0__2"
    // InternalElkGraph.g:1774:1: rule__ShapeLayout__Group_2_0__2 : rule__ShapeLayout__Group_2_0__2__Impl rule__ShapeLayout__Group_2_0__3 ;
    public final void rule__ShapeLayout__Group_2_0__2() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:1778:1: ( rule__ShapeLayout__Group_2_0__2__Impl rule__ShapeLayout__Group_2_0__3 )
            // InternalElkGraph.g:1779:2: rule__ShapeLayout__Group_2_0__2__Impl rule__ShapeLayout__Group_2_0__3
            {
            pushFollow(FOLLOW_18);
            rule__ShapeLayout__Group_2_0__2__Impl();

            state._fsp--;

            pushFollow(FOLLOW_2);
            rule__ShapeLayout__Group_2_0__3();

            state._fsp--;


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ShapeLayout__Group_2_0__2"


    // $ANTLR start "rule__ShapeLayout__Group_2_0__2__Impl"
    // InternalElkGraph.g:1786:1: rule__ShapeLayout__Group_2_0__2__Impl : ( ( rule__ShapeLayout__XAssignment_2_0_2 ) ) ;
    public final void rule__ShapeLayout__Group_2_0__2__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:1790:1: ( ( ( rule__ShapeLayout__XAssignment_2_0_2 ) ) )
            // InternalElkGraph.g:1791:1: ( ( rule__ShapeLayout__XAssignment_2_0_2 ) )
            {
            // InternalElkGraph.g:1791:1: ( ( rule__ShapeLayout__XAssignment_2_0_2 ) )
            // InternalElkGraph.g:1792:2: ( rule__ShapeLayout__XAssignment_2_0_2 )
            {
             before(grammarAccess.getShapeLayoutAccess().getXAssignment_2_0_2()); 
            // InternalElkGraph.g:1793:2: ( rule__ShapeLayout__XAssignment_2_0_2 )
            // InternalElkGraph.g:1793:3: rule__ShapeLayout__XAssignment_2_0_2
            {
            pushFollow(FOLLOW_2);
            rule__ShapeLayout__XAssignment_2_0_2();

            state._fsp--;


            }

             after(grammarAccess.getShapeLayoutAccess().getXAssignment_2_0_2()); 

            }


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ShapeLayout__Group_2_0__2__Impl"


    // $ANTLR start "rule__ShapeLayout__Group_2_0__3"
    // InternalElkGraph.g:1801:1: rule__ShapeLayout__Group_2_0__3 : rule__ShapeLayout__Group_2_0__3__Impl rule__ShapeLayout__Group_2_0__4 ;
    public final void rule__ShapeLayout__Group_2_0__3() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:1805:1: ( rule__ShapeLayout__Group_2_0__3__Impl rule__ShapeLayout__Group_2_0__4 )
            // InternalElkGraph.g:1806:2: rule__ShapeLayout__Group_2_0__3__Impl rule__ShapeLayout__Group_2_0__4
            {
            pushFollow(FOLLOW_17);
            rule__ShapeLayout__Group_2_0__3__Impl();

            state._fsp--;

            pushFollow(FOLLOW_2);
            rule__ShapeLayout__Group_2_0__4();

            state._fsp--;


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ShapeLayout__Group_2_0__3"


    // $ANTLR start "rule__ShapeLayout__Group_2_0__3__Impl"
    // InternalElkGraph.g:1813:1: rule__ShapeLayout__Group_2_0__3__Impl : ( ',' ) ;
    public final void rule__ShapeLayout__Group_2_0__3__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:1817:1: ( ( ',' ) )
            // InternalElkGraph.g:1818:1: ( ',' )
            {
            // InternalElkGraph.g:1818:1: ( ',' )
            // InternalElkGraph.g:1819:2: ','
            {
             before(grammarAccess.getShapeLayoutAccess().getCommaKeyword_2_0_3()); 
            match(input,27,FOLLOW_2); 
             after(grammarAccess.getShapeLayoutAccess().getCommaKeyword_2_0_3()); 

            }


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ShapeLayout__Group_2_0__3__Impl"


    // $ANTLR start "rule__ShapeLayout__Group_2_0__4"
    // InternalElkGraph.g:1828:1: rule__ShapeLayout__Group_2_0__4 : rule__ShapeLayout__Group_2_0__4__Impl ;
    public final void rule__ShapeLayout__Group_2_0__4() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:1832:1: ( rule__ShapeLayout__Group_2_0__4__Impl )
            // InternalElkGraph.g:1833:2: rule__ShapeLayout__Group_2_0__4__Impl
            {
            pushFollow(FOLLOW_2);
            rule__ShapeLayout__Group_2_0__4__Impl();

            state._fsp--;


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ShapeLayout__Group_2_0__4"


    // $ANTLR start "rule__ShapeLayout__Group_2_0__4__Impl"
    // InternalElkGraph.g:1839:1: rule__ShapeLayout__Group_2_0__4__Impl : ( ( rule__ShapeLayout__YAssignment_2_0_4 ) ) ;
    public final void rule__ShapeLayout__Group_2_0__4__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:1843:1: ( ( ( rule__ShapeLayout__YAssignment_2_0_4 ) ) )
            // InternalElkGraph.g:1844:1: ( ( rule__ShapeLayout__YAssignment_2_0_4 ) )
            {
            // InternalElkGraph.g:1844:1: ( ( rule__ShapeLayout__YAssignment_2_0_4 ) )
            // InternalElkGraph.g:1845:2: ( rule__ShapeLayout__YAssignment_2_0_4 )
            {
             before(grammarAccess.getShapeLayoutAccess().getYAssignment_2_0_4()); 
            // InternalElkGraph.g:1846:2: ( rule__ShapeLayout__YAssignment_2_0_4 )
            // InternalElkGraph.g:1846:3: rule__ShapeLayout__YAssignment_2_0_4
            {
            pushFollow(FOLLOW_2);
            rule__ShapeLayout__YAssignment_2_0_4();

            state._fsp--;


            }

             after(grammarAccess.getShapeLayoutAccess().getYAssignment_2_0_4()); 

            }


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ShapeLayout__Group_2_0__4__Impl"


    // $ANTLR start "rule__ShapeLayout__Group_2_1__0"
    // InternalElkGraph.g:1855:1: rule__ShapeLayout__Group_2_1__0 : rule__ShapeLayout__Group_2_1__0__Impl rule__ShapeLayout__Group_2_1__1 ;
    public final void rule__ShapeLayout__Group_2_1__0() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:1859:1: ( rule__ShapeLayout__Group_2_1__0__Impl rule__ShapeLayout__Group_2_1__1 )
            // InternalElkGraph.g:1860:2: rule__ShapeLayout__Group_2_1__0__Impl rule__ShapeLayout__Group_2_1__1
            {
            pushFollow(FOLLOW_11);
            rule__ShapeLayout__Group_2_1__0__Impl();

            state._fsp--;

            pushFollow(FOLLOW_2);
            rule__ShapeLayout__Group_2_1__1();

            state._fsp--;


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ShapeLayout__Group_2_1__0"


    // $ANTLR start "rule__ShapeLayout__Group_2_1__0__Impl"
    // InternalElkGraph.g:1867:1: rule__ShapeLayout__Group_2_1__0__Impl : ( 'size' ) ;
    public final void rule__ShapeLayout__Group_2_1__0__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:1871:1: ( ( 'size' ) )
            // InternalElkGraph.g:1872:1: ( 'size' )
            {
            // InternalElkGraph.g:1872:1: ( 'size' )
            // InternalElkGraph.g:1873:2: 'size'
            {
             before(grammarAccess.getShapeLayoutAccess().getSizeKeyword_2_1_0()); 
            match(input,28,FOLLOW_2); 
             after(grammarAccess.getShapeLayoutAccess().getSizeKeyword_2_1_0()); 

            }


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ShapeLayout__Group_2_1__0__Impl"


    // $ANTLR start "rule__ShapeLayout__Group_2_1__1"
    // InternalElkGraph.g:1882:1: rule__ShapeLayout__Group_2_1__1 : rule__ShapeLayout__Group_2_1__1__Impl rule__ShapeLayout__Group_2_1__2 ;
    public final void rule__ShapeLayout__Group_2_1__1() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:1886:1: ( rule__ShapeLayout__Group_2_1__1__Impl rule__ShapeLayout__Group_2_1__2 )
            // InternalElkGraph.g:1887:2: rule__ShapeLayout__Group_2_1__1__Impl rule__ShapeLayout__Group_2_1__2
            {
            pushFollow(FOLLOW_17);
            rule__ShapeLayout__Group_2_1__1__Impl();

            state._fsp--;

            pushFollow(FOLLOW_2);
            rule__ShapeLayout__Group_2_1__2();

            state._fsp--;


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ShapeLayout__Group_2_1__1"


    // $ANTLR start "rule__ShapeLayout__Group_2_1__1__Impl"
    // InternalElkGraph.g:1894:1: rule__ShapeLayout__Group_2_1__1__Impl : ( ':' ) ;
    public final void rule__ShapeLayout__Group_2_1__1__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:1898:1: ( ( ':' ) )
            // InternalElkGraph.g:1899:1: ( ':' )
            {
            // InternalElkGraph.g:1899:1: ( ':' )
            // InternalElkGraph.g:1900:2: ':'
            {
             before(grammarAccess.getShapeLayoutAccess().getColonKeyword_2_1_1()); 
            match(input,21,FOLLOW_2); 
             after(grammarAccess.getShapeLayoutAccess().getColonKeyword_2_1_1()); 

            }


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ShapeLayout__Group_2_1__1__Impl"


    // $ANTLR start "rule__ShapeLayout__Group_2_1__2"
    // InternalElkGraph.g:1909:1: rule__ShapeLayout__Group_2_1__2 : rule__ShapeLayout__Group_2_1__2__Impl rule__ShapeLayout__Group_2_1__3 ;
    public final void rule__ShapeLayout__Group_2_1__2() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:1913:1: ( rule__ShapeLayout__Group_2_1__2__Impl rule__ShapeLayout__Group_2_1__3 )
            // InternalElkGraph.g:1914:2: rule__ShapeLayout__Group_2_1__2__Impl rule__ShapeLayout__Group_2_1__3
            {
            pushFollow(FOLLOW_18);
            rule__ShapeLayout__Group_2_1__2__Impl();

            state._fsp--;

            pushFollow(FOLLOW_2);
            rule__ShapeLayout__Group_2_1__3();

            state._fsp--;


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ShapeLayout__Group_2_1__2"


    // $ANTLR start "rule__ShapeLayout__Group_2_1__2__Impl"
    // InternalElkGraph.g:1921:1: rule__ShapeLayout__Group_2_1__2__Impl : ( ( rule__ShapeLayout__WidthAssignment_2_1_2 ) ) ;
    public final void rule__ShapeLayout__Group_2_1__2__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:1925:1: ( ( ( rule__ShapeLayout__WidthAssignment_2_1_2 ) ) )
            // InternalElkGraph.g:1926:1: ( ( rule__ShapeLayout__WidthAssignment_2_1_2 ) )
            {
            // InternalElkGraph.g:1926:1: ( ( rule__ShapeLayout__WidthAssignment_2_1_2 ) )
            // InternalElkGraph.g:1927:2: ( rule__ShapeLayout__WidthAssignment_2_1_2 )
            {
             before(grammarAccess.getShapeLayoutAccess().getWidthAssignment_2_1_2()); 
            // InternalElkGraph.g:1928:2: ( rule__ShapeLayout__WidthAssignment_2_1_2 )
            // InternalElkGraph.g:1928:3: rule__ShapeLayout__WidthAssignment_2_1_2
            {
            pushFollow(FOLLOW_2);
            rule__ShapeLayout__WidthAssignment_2_1_2();

            state._fsp--;


            }

             after(grammarAccess.getShapeLayoutAccess().getWidthAssignment_2_1_2()); 

            }


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ShapeLayout__Group_2_1__2__Impl"


    // $ANTLR start "rule__ShapeLayout__Group_2_1__3"
    // InternalElkGraph.g:1936:1: rule__ShapeLayout__Group_2_1__3 : rule__ShapeLayout__Group_2_1__3__Impl rule__ShapeLayout__Group_2_1__4 ;
    public final void rule__ShapeLayout__Group_2_1__3() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:1940:1: ( rule__ShapeLayout__Group_2_1__3__Impl rule__ShapeLayout__Group_2_1__4 )
            // InternalElkGraph.g:1941:2: rule__ShapeLayout__Group_2_1__3__Impl rule__ShapeLayout__Group_2_1__4
            {
            pushFollow(FOLLOW_17);
            rule__ShapeLayout__Group_2_1__3__Impl();

            state._fsp--;

            pushFollow(FOLLOW_2);
            rule__ShapeLayout__Group_2_1__4();

            state._fsp--;


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ShapeLayout__Group_2_1__3"


    // $ANTLR start "rule__ShapeLayout__Group_2_1__3__Impl"
    // InternalElkGraph.g:1948:1: rule__ShapeLayout__Group_2_1__3__Impl : ( ',' ) ;
    public final void rule__ShapeLayout__Group_2_1__3__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:1952:1: ( ( ',' ) )
            // InternalElkGraph.g:1953:1: ( ',' )
            {
            // InternalElkGraph.g:1953:1: ( ',' )
            // InternalElkGraph.g:1954:2: ','
            {
             before(grammarAccess.getShapeLayoutAccess().getCommaKeyword_2_1_3()); 
            match(input,27,FOLLOW_2); 
             after(grammarAccess.getShapeLayoutAccess().getCommaKeyword_2_1_3()); 

            }


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ShapeLayout__Group_2_1__3__Impl"


    // $ANTLR start "rule__ShapeLayout__Group_2_1__4"
    // InternalElkGraph.g:1963:1: rule__ShapeLayout__Group_2_1__4 : rule__ShapeLayout__Group_2_1__4__Impl ;
    public final void rule__ShapeLayout__Group_2_1__4() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:1967:1: ( rule__ShapeLayout__Group_2_1__4__Impl )
            // InternalElkGraph.g:1968:2: rule__ShapeLayout__Group_2_1__4__Impl
            {
            pushFollow(FOLLOW_2);
            rule__ShapeLayout__Group_2_1__4__Impl();

            state._fsp--;


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ShapeLayout__Group_2_1__4"


    // $ANTLR start "rule__ShapeLayout__Group_2_1__4__Impl"
    // InternalElkGraph.g:1974:1: rule__ShapeLayout__Group_2_1__4__Impl : ( ( rule__ShapeLayout__HeightAssignment_2_1_4 ) ) ;
    public final void rule__ShapeLayout__Group_2_1__4__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:1978:1: ( ( ( rule__ShapeLayout__HeightAssignment_2_1_4 ) ) )
            // InternalElkGraph.g:1979:1: ( ( rule__ShapeLayout__HeightAssignment_2_1_4 ) )
            {
            // InternalElkGraph.g:1979:1: ( ( rule__ShapeLayout__HeightAssignment_2_1_4 ) )
            // InternalElkGraph.g:1980:2: ( rule__ShapeLayout__HeightAssignment_2_1_4 )
            {
             before(grammarAccess.getShapeLayoutAccess().getHeightAssignment_2_1_4()); 
            // InternalElkGraph.g:1981:2: ( rule__ShapeLayout__HeightAssignment_2_1_4 )
            // InternalElkGraph.g:1981:3: rule__ShapeLayout__HeightAssignment_2_1_4
            {
            pushFollow(FOLLOW_2);
            rule__ShapeLayout__HeightAssignment_2_1_4();

            state._fsp--;


            }

             after(grammarAccess.getShapeLayoutAccess().getHeightAssignment_2_1_4()); 

            }


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ShapeLayout__Group_2_1__4__Impl"


    // $ANTLR start "rule__ElkEdge__Group__0"
    // InternalElkGraph.g:1990:1: rule__ElkEdge__Group__0 : rule__ElkEdge__Group__0__Impl rule__ElkEdge__Group__1 ;
    public final void rule__ElkEdge__Group__0() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:1994:1: ( rule__ElkEdge__Group__0__Impl rule__ElkEdge__Group__1 )
            // InternalElkGraph.g:1995:2: rule__ElkEdge__Group__0__Impl rule__ElkEdge__Group__1
            {
            pushFollow(FOLLOW_7);
            rule__ElkEdge__Group__0__Impl();

            state._fsp--;

            pushFollow(FOLLOW_2);
            rule__ElkEdge__Group__1();

            state._fsp--;


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ElkEdge__Group__0"


    // $ANTLR start "rule__ElkEdge__Group__0__Impl"
    // InternalElkGraph.g:2002:1: rule__ElkEdge__Group__0__Impl : ( 'edge' ) ;
    public final void rule__ElkEdge__Group__0__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:2006:1: ( ( 'edge' ) )
            // InternalElkGraph.g:2007:1: ( 'edge' )
            {
            // InternalElkGraph.g:2007:1: ( 'edge' )
            // InternalElkGraph.g:2008:2: 'edge'
            {
             before(grammarAccess.getElkEdgeAccess().getEdgeKeyword_0()); 
            match(input,29,FOLLOW_2); 
             after(grammarAccess.getElkEdgeAccess().getEdgeKeyword_0()); 

            }


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ElkEdge__Group__0__Impl"


    // $ANTLR start "rule__ElkEdge__Group__1"
    // InternalElkGraph.g:2017:1: rule__ElkEdge__Group__1 : rule__ElkEdge__Group__1__Impl rule__ElkEdge__Group__2 ;
    public final void rule__ElkEdge__Group__1() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:2021:1: ( rule__ElkEdge__Group__1__Impl rule__ElkEdge__Group__2 )
            // InternalElkGraph.g:2022:2: rule__ElkEdge__Group__1__Impl rule__ElkEdge__Group__2
            {
            pushFollow(FOLLOW_7);
            rule__ElkEdge__Group__1__Impl();

            state._fsp--;

            pushFollow(FOLLOW_2);
            rule__ElkEdge__Group__2();

            state._fsp--;


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ElkEdge__Group__1"


    // $ANTLR start "rule__ElkEdge__Group__1__Impl"
    // InternalElkGraph.g:2029:1: rule__ElkEdge__Group__1__Impl : ( ( rule__ElkEdge__Group_1__0 )? ) ;
    public final void rule__ElkEdge__Group__1__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:2033:1: ( ( ( rule__ElkEdge__Group_1__0 )? ) )
            // InternalElkGraph.g:2034:1: ( ( rule__ElkEdge__Group_1__0 )? )
            {
            // InternalElkGraph.g:2034:1: ( ( rule__ElkEdge__Group_1__0 )? )
            // InternalElkGraph.g:2035:2: ( rule__ElkEdge__Group_1__0 )?
            {
             before(grammarAccess.getElkEdgeAccess().getGroup_1()); 
            // InternalElkGraph.g:2036:2: ( rule__ElkEdge__Group_1__0 )?
            int alt26=2;
            int LA26_0 = input.LA(1);

            if ( (LA26_0==RULE_ID) ) {
                int LA26_1 = input.LA(2);

                if ( (LA26_1==21) ) {
                    alt26=1;
                }
            }
            switch (alt26) {
                case 1 :
                    // InternalElkGraph.g:2036:3: rule__ElkEdge__Group_1__0
                    {
                    pushFollow(FOLLOW_2);
                    rule__ElkEdge__Group_1__0();

                    state._fsp--;


                    }
                    break;

            }

             after(grammarAccess.getElkEdgeAccess().getGroup_1()); 

            }


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ElkEdge__Group__1__Impl"


    // $ANTLR start "rule__ElkEdge__Group__2"
    // InternalElkGraph.g:2044:1: rule__ElkEdge__Group__2 : rule__ElkEdge__Group__2__Impl rule__ElkEdge__Group__3 ;
    public final void rule__ElkEdge__Group__2() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:2048:1: ( rule__ElkEdge__Group__2__Impl rule__ElkEdge__Group__3 )
            // InternalElkGraph.g:2049:2: rule__ElkEdge__Group__2__Impl rule__ElkEdge__Group__3
            {
            pushFollow(FOLLOW_19);
            rule__ElkEdge__Group__2__Impl();

            state._fsp--;

            pushFollow(FOLLOW_2);
            rule__ElkEdge__Group__3();

            state._fsp--;


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ElkEdge__Group__2"


    // $ANTLR start "rule__ElkEdge__Group__2__Impl"
    // InternalElkGraph.g:2056:1: rule__ElkEdge__Group__2__Impl : ( ( rule__ElkEdge__SourcesAssignment_2 ) ) ;
    public final void rule__ElkEdge__Group__2__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:2060:1: ( ( ( rule__ElkEdge__SourcesAssignment_2 ) ) )
            // InternalElkGraph.g:2061:1: ( ( rule__ElkEdge__SourcesAssignment_2 ) )
            {
            // InternalElkGraph.g:2061:1: ( ( rule__ElkEdge__SourcesAssignment_2 ) )
            // InternalElkGraph.g:2062:2: ( rule__ElkEdge__SourcesAssignment_2 )
            {
             before(grammarAccess.getElkEdgeAccess().getSourcesAssignment_2()); 
            // InternalElkGraph.g:2063:2: ( rule__ElkEdge__SourcesAssignment_2 )
            // InternalElkGraph.g:2063:3: rule__ElkEdge__SourcesAssignment_2
            {
            pushFollow(FOLLOW_2);
            rule__ElkEdge__SourcesAssignment_2();

            state._fsp--;


            }

             after(grammarAccess.getElkEdgeAccess().getSourcesAssignment_2()); 

            }


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ElkEdge__Group__2__Impl"


    // $ANTLR start "rule__ElkEdge__Group__3"
    // InternalElkGraph.g:2071:1: rule__ElkEdge__Group__3 : rule__ElkEdge__Group__3__Impl rule__ElkEdge__Group__4 ;
    public final void rule__ElkEdge__Group__3() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:2075:1: ( rule__ElkEdge__Group__3__Impl rule__ElkEdge__Group__4 )
            // InternalElkGraph.g:2076:2: rule__ElkEdge__Group__3__Impl rule__ElkEdge__Group__4
            {
            pushFollow(FOLLOW_19);
            rule__ElkEdge__Group__3__Impl();

            state._fsp--;

            pushFollow(FOLLOW_2);
            rule__ElkEdge__Group__4();

            state._fsp--;


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ElkEdge__Group__3"


    // $ANTLR start "rule__ElkEdge__Group__3__Impl"
    // InternalElkGraph.g:2083:1: rule__ElkEdge__Group__3__Impl : ( ( rule__ElkEdge__Group_3__0 )* ) ;
    public final void rule__ElkEdge__Group__3__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:2087:1: ( ( ( rule__ElkEdge__Group_3__0 )* ) )
            // InternalElkGraph.g:2088:1: ( ( rule__ElkEdge__Group_3__0 )* )
            {
            // InternalElkGraph.g:2088:1: ( ( rule__ElkEdge__Group_3__0 )* )
            // InternalElkGraph.g:2089:2: ( rule__ElkEdge__Group_3__0 )*
            {
             before(grammarAccess.getElkEdgeAccess().getGroup_3()); 
            // InternalElkGraph.g:2090:2: ( rule__ElkEdge__Group_3__0 )*
            loop27:
            do {
                int alt27=2;
                int LA27_0 = input.LA(1);

                if ( (LA27_0==27) ) {
                    alt27=1;
                }


                switch (alt27) {
            	case 1 :
            	    // InternalElkGraph.g:2090:3: rule__ElkEdge__Group_3__0
            	    {
            	    pushFollow(FOLLOW_20);
            	    rule__ElkEdge__Group_3__0();

            	    state._fsp--;


            	    }
            	    break;

            	default :
            	    break loop27;
                }
            } while (true);

             after(grammarAccess.getElkEdgeAccess().getGroup_3()); 

            }


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ElkEdge__Group__3__Impl"


    // $ANTLR start "rule__ElkEdge__Group__4"
    // InternalElkGraph.g:2098:1: rule__ElkEdge__Group__4 : rule__ElkEdge__Group__4__Impl rule__ElkEdge__Group__5 ;
    public final void rule__ElkEdge__Group__4() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:2102:1: ( rule__ElkEdge__Group__4__Impl rule__ElkEdge__Group__5 )
            // InternalElkGraph.g:2103:2: rule__ElkEdge__Group__4__Impl rule__ElkEdge__Group__5
            {
            pushFollow(FOLLOW_7);
            rule__ElkEdge__Group__4__Impl();

            state._fsp--;

            pushFollow(FOLLOW_2);
            rule__ElkEdge__Group__5();

            state._fsp--;


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ElkEdge__Group__4"


    // $ANTLR start "rule__ElkEdge__Group__4__Impl"
    // InternalElkGraph.g:2110:1: rule__ElkEdge__Group__4__Impl : ( '->' ) ;
    public final void rule__ElkEdge__Group__4__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:2114:1: ( ( '->' ) )
            // InternalElkGraph.g:2115:1: ( '->' )
            {
            // InternalElkGraph.g:2115:1: ( '->' )
            // InternalElkGraph.g:2116:2: '->'
            {
             before(grammarAccess.getElkEdgeAccess().getHyphenMinusGreaterThanSignKeyword_4()); 
            match(input,30,FOLLOW_2); 
             after(grammarAccess.getElkEdgeAccess().getHyphenMinusGreaterThanSignKeyword_4()); 

            }


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ElkEdge__Group__4__Impl"


    // $ANTLR start "rule__ElkEdge__Group__5"
    // InternalElkGraph.g:2125:1: rule__ElkEdge__Group__5 : rule__ElkEdge__Group__5__Impl rule__ElkEdge__Group__6 ;
    public final void rule__ElkEdge__Group__5() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:2129:1: ( rule__ElkEdge__Group__5__Impl rule__ElkEdge__Group__6 )
            // InternalElkGraph.g:2130:2: rule__ElkEdge__Group__5__Impl rule__ElkEdge__Group__6
            {
            pushFollow(FOLLOW_21);
            rule__ElkEdge__Group__5__Impl();

            state._fsp--;

            pushFollow(FOLLOW_2);
            rule__ElkEdge__Group__6();

            state._fsp--;


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ElkEdge__Group__5"


    // $ANTLR start "rule__ElkEdge__Group__5__Impl"
    // InternalElkGraph.g:2137:1: rule__ElkEdge__Group__5__Impl : ( ( rule__ElkEdge__TargetsAssignment_5 ) ) ;
    public final void rule__ElkEdge__Group__5__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:2141:1: ( ( ( rule__ElkEdge__TargetsAssignment_5 ) ) )
            // InternalElkGraph.g:2142:1: ( ( rule__ElkEdge__TargetsAssignment_5 ) )
            {
            // InternalElkGraph.g:2142:1: ( ( rule__ElkEdge__TargetsAssignment_5 ) )
            // InternalElkGraph.g:2143:2: ( rule__ElkEdge__TargetsAssignment_5 )
            {
             before(grammarAccess.getElkEdgeAccess().getTargetsAssignment_5()); 
            // InternalElkGraph.g:2144:2: ( rule__ElkEdge__TargetsAssignment_5 )
            // InternalElkGraph.g:2144:3: rule__ElkEdge__TargetsAssignment_5
            {
            pushFollow(FOLLOW_2);
            rule__ElkEdge__TargetsAssignment_5();

            state._fsp--;


            }

             after(grammarAccess.getElkEdgeAccess().getTargetsAssignment_5()); 

            }


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ElkEdge__Group__5__Impl"


    // $ANTLR start "rule__ElkEdge__Group__6"
    // InternalElkGraph.g:2152:1: rule__ElkEdge__Group__6 : rule__ElkEdge__Group__6__Impl rule__ElkEdge__Group__7 ;
    public final void rule__ElkEdge__Group__6() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:2156:1: ( rule__ElkEdge__Group__6__Impl rule__ElkEdge__Group__7 )
            // InternalElkGraph.g:2157:2: rule__ElkEdge__Group__6__Impl rule__ElkEdge__Group__7
            {
            pushFollow(FOLLOW_21);
            rule__ElkEdge__Group__6__Impl();

            state._fsp--;

            pushFollow(FOLLOW_2);
            rule__ElkEdge__Group__7();

            state._fsp--;


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ElkEdge__Group__6"


    // $ANTLR start "rule__ElkEdge__Group__6__Impl"
    // InternalElkGraph.g:2164:1: rule__ElkEdge__Group__6__Impl : ( ( rule__ElkEdge__Group_6__0 )* ) ;
    public final void rule__ElkEdge__Group__6__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:2168:1: ( ( ( rule__ElkEdge__Group_6__0 )* ) )
            // InternalElkGraph.g:2169:1: ( ( rule__ElkEdge__Group_6__0 )* )
            {
            // InternalElkGraph.g:2169:1: ( ( rule__ElkEdge__Group_6__0 )* )
            // InternalElkGraph.g:2170:2: ( rule__ElkEdge__Group_6__0 )*
            {
             before(grammarAccess.getElkEdgeAccess().getGroup_6()); 
            // InternalElkGraph.g:2171:2: ( rule__ElkEdge__Group_6__0 )*
            loop28:
            do {
                int alt28=2;
                int LA28_0 = input.LA(1);

                if ( (LA28_0==27) ) {
                    alt28=1;
                }


                switch (alt28) {
            	case 1 :
            	    // InternalElkGraph.g:2171:3: rule__ElkEdge__Group_6__0
            	    {
            	    pushFollow(FOLLOW_20);
            	    rule__ElkEdge__Group_6__0();

            	    state._fsp--;


            	    }
            	    break;

            	default :
            	    break loop28;
                }
            } while (true);

             after(grammarAccess.getElkEdgeAccess().getGroup_6()); 

            }


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ElkEdge__Group__6__Impl"


    // $ANTLR start "rule__ElkEdge__Group__7"
    // InternalElkGraph.g:2179:1: rule__ElkEdge__Group__7 : rule__ElkEdge__Group__7__Impl ;
    public final void rule__ElkEdge__Group__7() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:2183:1: ( rule__ElkEdge__Group__7__Impl )
            // InternalElkGraph.g:2184:2: rule__ElkEdge__Group__7__Impl
            {
            pushFollow(FOLLOW_2);
            rule__ElkEdge__Group__7__Impl();

            state._fsp--;


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ElkEdge__Group__7"


    // $ANTLR start "rule__ElkEdge__Group__7__Impl"
    // InternalElkGraph.g:2190:1: rule__ElkEdge__Group__7__Impl : ( ( rule__ElkEdge__Group_7__0 )? ) ;
    public final void rule__ElkEdge__Group__7__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:2194:1: ( ( ( rule__ElkEdge__Group_7__0 )? ) )
            // InternalElkGraph.g:2195:1: ( ( rule__ElkEdge__Group_7__0 )? )
            {
            // InternalElkGraph.g:2195:1: ( ( rule__ElkEdge__Group_7__0 )? )
            // InternalElkGraph.g:2196:2: ( rule__ElkEdge__Group_7__0 )?
            {
             before(grammarAccess.getElkEdgeAccess().getGroup_7()); 
            // InternalElkGraph.g:2197:2: ( rule__ElkEdge__Group_7__0 )?
            int alt29=2;
            int LA29_0 = input.LA(1);

            if ( (LA29_0==18) ) {
                alt29=1;
            }
            switch (alt29) {
                case 1 :
                    // InternalElkGraph.g:2197:3: rule__ElkEdge__Group_7__0
                    {
                    pushFollow(FOLLOW_2);
                    rule__ElkEdge__Group_7__0();

                    state._fsp--;


                    }
                    break;

            }

             after(grammarAccess.getElkEdgeAccess().getGroup_7()); 

            }


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ElkEdge__Group__7__Impl"


    // $ANTLR start "rule__ElkEdge__Group_1__0"
    // InternalElkGraph.g:2206:1: rule__ElkEdge__Group_1__0 : rule__ElkEdge__Group_1__0__Impl rule__ElkEdge__Group_1__1 ;
    public final void rule__ElkEdge__Group_1__0() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:2210:1: ( rule__ElkEdge__Group_1__0__Impl rule__ElkEdge__Group_1__1 )
            // InternalElkGraph.g:2211:2: rule__ElkEdge__Group_1__0__Impl rule__ElkEdge__Group_1__1
            {
            pushFollow(FOLLOW_11);
            rule__ElkEdge__Group_1__0__Impl();

            state._fsp--;

            pushFollow(FOLLOW_2);
            rule__ElkEdge__Group_1__1();

            state._fsp--;


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ElkEdge__Group_1__0"


    // $ANTLR start "rule__ElkEdge__Group_1__0__Impl"
    // InternalElkGraph.g:2218:1: rule__ElkEdge__Group_1__0__Impl : ( ( rule__ElkEdge__IdentifierAssignment_1_0 ) ) ;
    public final void rule__ElkEdge__Group_1__0__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:2222:1: ( ( ( rule__ElkEdge__IdentifierAssignment_1_0 ) ) )
            // InternalElkGraph.g:2223:1: ( ( rule__ElkEdge__IdentifierAssignment_1_0 ) )
            {
            // InternalElkGraph.g:2223:1: ( ( rule__ElkEdge__IdentifierAssignment_1_0 ) )
            // InternalElkGraph.g:2224:2: ( rule__ElkEdge__IdentifierAssignment_1_0 )
            {
             before(grammarAccess.getElkEdgeAccess().getIdentifierAssignment_1_0()); 
            // InternalElkGraph.g:2225:2: ( rule__ElkEdge__IdentifierAssignment_1_0 )
            // InternalElkGraph.g:2225:3: rule__ElkEdge__IdentifierAssignment_1_0
            {
            pushFollow(FOLLOW_2);
            rule__ElkEdge__IdentifierAssignment_1_0();

            state._fsp--;


            }

             after(grammarAccess.getElkEdgeAccess().getIdentifierAssignment_1_0()); 

            }


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ElkEdge__Group_1__0__Impl"


    // $ANTLR start "rule__ElkEdge__Group_1__1"
    // InternalElkGraph.g:2233:1: rule__ElkEdge__Group_1__1 : rule__ElkEdge__Group_1__1__Impl ;
    public final void rule__ElkEdge__Group_1__1() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:2237:1: ( rule__ElkEdge__Group_1__1__Impl )
            // InternalElkGraph.g:2238:2: rule__ElkEdge__Group_1__1__Impl
            {
            pushFollow(FOLLOW_2);
            rule__ElkEdge__Group_1__1__Impl();

            state._fsp--;


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ElkEdge__Group_1__1"


    // $ANTLR start "rule__ElkEdge__Group_1__1__Impl"
    // InternalElkGraph.g:2244:1: rule__ElkEdge__Group_1__1__Impl : ( ':' ) ;
    public final void rule__ElkEdge__Group_1__1__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:2248:1: ( ( ':' ) )
            // InternalElkGraph.g:2249:1: ( ':' )
            {
            // InternalElkGraph.g:2249:1: ( ':' )
            // InternalElkGraph.g:2250:2: ':'
            {
             before(grammarAccess.getElkEdgeAccess().getColonKeyword_1_1()); 
            match(input,21,FOLLOW_2); 
             after(grammarAccess.getElkEdgeAccess().getColonKeyword_1_1()); 

            }


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ElkEdge__Group_1__1__Impl"


    // $ANTLR start "rule__ElkEdge__Group_3__0"
    // InternalElkGraph.g:2260:1: rule__ElkEdge__Group_3__0 : rule__ElkEdge__Group_3__0__Impl rule__ElkEdge__Group_3__1 ;
    public final void rule__ElkEdge__Group_3__0() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:2264:1: ( rule__ElkEdge__Group_3__0__Impl rule__ElkEdge__Group_3__1 )
            // InternalElkGraph.g:2265:2: rule__ElkEdge__Group_3__0__Impl rule__ElkEdge__Group_3__1
            {
            pushFollow(FOLLOW_7);
            rule__ElkEdge__Group_3__0__Impl();

            state._fsp--;

            pushFollow(FOLLOW_2);
            rule__ElkEdge__Group_3__1();

            state._fsp--;


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ElkEdge__Group_3__0"


    // $ANTLR start "rule__ElkEdge__Group_3__0__Impl"
    // InternalElkGraph.g:2272:1: rule__ElkEdge__Group_3__0__Impl : ( ',' ) ;
    public final void rule__ElkEdge__Group_3__0__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:2276:1: ( ( ',' ) )
            // InternalElkGraph.g:2277:1: ( ',' )
            {
            // InternalElkGraph.g:2277:1: ( ',' )
            // InternalElkGraph.g:2278:2: ','
            {
             before(grammarAccess.getElkEdgeAccess().getCommaKeyword_3_0()); 
            match(input,27,FOLLOW_2); 
             after(grammarAccess.getElkEdgeAccess().getCommaKeyword_3_0()); 

            }


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ElkEdge__Group_3__0__Impl"


    // $ANTLR start "rule__ElkEdge__Group_3__1"
    // InternalElkGraph.g:2287:1: rule__ElkEdge__Group_3__1 : rule__ElkEdge__Group_3__1__Impl ;
    public final void rule__ElkEdge__Group_3__1() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:2291:1: ( rule__ElkEdge__Group_3__1__Impl )
            // InternalElkGraph.g:2292:2: rule__ElkEdge__Group_3__1__Impl
            {
            pushFollow(FOLLOW_2);
            rule__ElkEdge__Group_3__1__Impl();

            state._fsp--;


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ElkEdge__Group_3__1"


    // $ANTLR start "rule__ElkEdge__Group_3__1__Impl"
    // InternalElkGraph.g:2298:1: rule__ElkEdge__Group_3__1__Impl : ( ( rule__ElkEdge__SourcesAssignment_3_1 ) ) ;
    public final void rule__ElkEdge__Group_3__1__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:2302:1: ( ( ( rule__ElkEdge__SourcesAssignment_3_1 ) ) )
            // InternalElkGraph.g:2303:1: ( ( rule__ElkEdge__SourcesAssignment_3_1 ) )
            {
            // InternalElkGraph.g:2303:1: ( ( rule__ElkEdge__SourcesAssignment_3_1 ) )
            // InternalElkGraph.g:2304:2: ( rule__ElkEdge__SourcesAssignment_3_1 )
            {
             before(grammarAccess.getElkEdgeAccess().getSourcesAssignment_3_1()); 
            // InternalElkGraph.g:2305:2: ( rule__ElkEdge__SourcesAssignment_3_1 )
            // InternalElkGraph.g:2305:3: rule__ElkEdge__SourcesAssignment_3_1
            {
            pushFollow(FOLLOW_2);
            rule__ElkEdge__SourcesAssignment_3_1();

            state._fsp--;


            }

             after(grammarAccess.getElkEdgeAccess().getSourcesAssignment_3_1()); 

            }


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ElkEdge__Group_3__1__Impl"


    // $ANTLR start "rule__ElkEdge__Group_6__0"
    // InternalElkGraph.g:2314:1: rule__ElkEdge__Group_6__0 : rule__ElkEdge__Group_6__0__Impl rule__ElkEdge__Group_6__1 ;
    public final void rule__ElkEdge__Group_6__0() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:2318:1: ( rule__ElkEdge__Group_6__0__Impl rule__ElkEdge__Group_6__1 )
            // InternalElkGraph.g:2319:2: rule__ElkEdge__Group_6__0__Impl rule__ElkEdge__Group_6__1
            {
            pushFollow(FOLLOW_7);
            rule__ElkEdge__Group_6__0__Impl();

            state._fsp--;

            pushFollow(FOLLOW_2);
            rule__ElkEdge__Group_6__1();

            state._fsp--;


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ElkEdge__Group_6__0"


    // $ANTLR start "rule__ElkEdge__Group_6__0__Impl"
    // InternalElkGraph.g:2326:1: rule__ElkEdge__Group_6__0__Impl : ( ',' ) ;
    public final void rule__ElkEdge__Group_6__0__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:2330:1: ( ( ',' ) )
            // InternalElkGraph.g:2331:1: ( ',' )
            {
            // InternalElkGraph.g:2331:1: ( ',' )
            // InternalElkGraph.g:2332:2: ','
            {
             before(grammarAccess.getElkEdgeAccess().getCommaKeyword_6_0()); 
            match(input,27,FOLLOW_2); 
             after(grammarAccess.getElkEdgeAccess().getCommaKeyword_6_0()); 

            }


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ElkEdge__Group_6__0__Impl"


    // $ANTLR start "rule__ElkEdge__Group_6__1"
    // InternalElkGraph.g:2341:1: rule__ElkEdge__Group_6__1 : rule__ElkEdge__Group_6__1__Impl ;
    public final void rule__ElkEdge__Group_6__1() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:2345:1: ( rule__ElkEdge__Group_6__1__Impl )
            // InternalElkGraph.g:2346:2: rule__ElkEdge__Group_6__1__Impl
            {
            pushFollow(FOLLOW_2);
            rule__ElkEdge__Group_6__1__Impl();

            state._fsp--;


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ElkEdge__Group_6__1"


    // $ANTLR start "rule__ElkEdge__Group_6__1__Impl"
    // InternalElkGraph.g:2352:1: rule__ElkEdge__Group_6__1__Impl : ( ( rule__ElkEdge__TargetsAssignment_6_1 ) ) ;
    public final void rule__ElkEdge__Group_6__1__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:2356:1: ( ( ( rule__ElkEdge__TargetsAssignment_6_1 ) ) )
            // InternalElkGraph.g:2357:1: ( ( rule__ElkEdge__TargetsAssignment_6_1 ) )
            {
            // InternalElkGraph.g:2357:1: ( ( rule__ElkEdge__TargetsAssignment_6_1 ) )
            // InternalElkGraph.g:2358:2: ( rule__ElkEdge__TargetsAssignment_6_1 )
            {
             before(grammarAccess.getElkEdgeAccess().getTargetsAssignment_6_1()); 
            // InternalElkGraph.g:2359:2: ( rule__ElkEdge__TargetsAssignment_6_1 )
            // InternalElkGraph.g:2359:3: rule__ElkEdge__TargetsAssignment_6_1
            {
            pushFollow(FOLLOW_2);
            rule__ElkEdge__TargetsAssignment_6_1();

            state._fsp--;


            }

             after(grammarAccess.getElkEdgeAccess().getTargetsAssignment_6_1()); 

            }


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ElkEdge__Group_6__1__Impl"


    // $ANTLR start "rule__ElkEdge__Group_7__0"
    // InternalElkGraph.g:2368:1: rule__ElkEdge__Group_7__0 : rule__ElkEdge__Group_7__0__Impl rule__ElkEdge__Group_7__1 ;
    public final void rule__ElkEdge__Group_7__0() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:2372:1: ( rule__ElkEdge__Group_7__0__Impl rule__ElkEdge__Group_7__1 )
            // InternalElkGraph.g:2373:2: rule__ElkEdge__Group_7__0__Impl rule__ElkEdge__Group_7__1
            {
            pushFollow(FOLLOW_12);
            rule__ElkEdge__Group_7__0__Impl();

            state._fsp--;

            pushFollow(FOLLOW_2);
            rule__ElkEdge__Group_7__1();

            state._fsp--;


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ElkEdge__Group_7__0"


    // $ANTLR start "rule__ElkEdge__Group_7__0__Impl"
    // InternalElkGraph.g:2380:1: rule__ElkEdge__Group_7__0__Impl : ( '{' ) ;
    public final void rule__ElkEdge__Group_7__0__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:2384:1: ( ( '{' ) )
            // InternalElkGraph.g:2385:1: ( '{' )
            {
            // InternalElkGraph.g:2385:1: ( '{' )
            // InternalElkGraph.g:2386:2: '{'
            {
             before(grammarAccess.getElkEdgeAccess().getLeftCurlyBracketKeyword_7_0()); 
            match(input,18,FOLLOW_2); 
             after(grammarAccess.getElkEdgeAccess().getLeftCurlyBracketKeyword_7_0()); 

            }


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ElkEdge__Group_7__0__Impl"


    // $ANTLR start "rule__ElkEdge__Group_7__1"
    // InternalElkGraph.g:2395:1: rule__ElkEdge__Group_7__1 : rule__ElkEdge__Group_7__1__Impl rule__ElkEdge__Group_7__2 ;
    public final void rule__ElkEdge__Group_7__1() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:2399:1: ( rule__ElkEdge__Group_7__1__Impl rule__ElkEdge__Group_7__2 )
            // InternalElkGraph.g:2400:2: rule__ElkEdge__Group_7__1__Impl rule__ElkEdge__Group_7__2
            {
            pushFollow(FOLLOW_12);
            rule__ElkEdge__Group_7__1__Impl();

            state._fsp--;

            pushFollow(FOLLOW_2);
            rule__ElkEdge__Group_7__2();

            state._fsp--;


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ElkEdge__Group_7__1"


    // $ANTLR start "rule__ElkEdge__Group_7__1__Impl"
    // InternalElkGraph.g:2407:1: rule__ElkEdge__Group_7__1__Impl : ( ( ruleEdgeLayout )? ) ;
    public final void rule__ElkEdge__Group_7__1__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:2411:1: ( ( ( ruleEdgeLayout )? ) )
            // InternalElkGraph.g:2412:1: ( ( ruleEdgeLayout )? )
            {
            // InternalElkGraph.g:2412:1: ( ( ruleEdgeLayout )? )
            // InternalElkGraph.g:2413:2: ( ruleEdgeLayout )?
            {
             before(grammarAccess.getElkEdgeAccess().getEdgeLayoutParserRuleCall_7_1()); 
            // InternalElkGraph.g:2414:2: ( ruleEdgeLayout )?
            int alt30=2;
            int LA30_0 = input.LA(1);

            if ( (LA30_0==23) ) {
                alt30=1;
            }
            switch (alt30) {
                case 1 :
                    // InternalElkGraph.g:2414:3: ruleEdgeLayout
                    {
                    pushFollow(FOLLOW_2);
                    ruleEdgeLayout();

                    state._fsp--;


                    }
                    break;

            }

             after(grammarAccess.getElkEdgeAccess().getEdgeLayoutParserRuleCall_7_1()); 

            }


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ElkEdge__Group_7__1__Impl"


    // $ANTLR start "rule__ElkEdge__Group_7__2"
    // InternalElkGraph.g:2422:1: rule__ElkEdge__Group_7__2 : rule__ElkEdge__Group_7__2__Impl rule__ElkEdge__Group_7__3 ;
    public final void rule__ElkEdge__Group_7__2() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:2426:1: ( rule__ElkEdge__Group_7__2__Impl rule__ElkEdge__Group_7__3 )
            // InternalElkGraph.g:2427:2: rule__ElkEdge__Group_7__2__Impl rule__ElkEdge__Group_7__3
            {
            pushFollow(FOLLOW_12);
            rule__ElkEdge__Group_7__2__Impl();

            state._fsp--;

            pushFollow(FOLLOW_2);
            rule__ElkEdge__Group_7__3();

            state._fsp--;


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ElkEdge__Group_7__2"


    // $ANTLR start "rule__ElkEdge__Group_7__2__Impl"
    // InternalElkGraph.g:2434:1: rule__ElkEdge__Group_7__2__Impl : ( ( rule__ElkEdge__PropertiesAssignment_7_2 )* ) ;
    public final void rule__ElkEdge__Group_7__2__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:2438:1: ( ( ( rule__ElkEdge__PropertiesAssignment_7_2 )* ) )
            // InternalElkGraph.g:2439:1: ( ( rule__ElkEdge__PropertiesAssignment_7_2 )* )
            {
            // InternalElkGraph.g:2439:1: ( ( rule__ElkEdge__PropertiesAssignment_7_2 )* )
            // InternalElkGraph.g:2440:2: ( rule__ElkEdge__PropertiesAssignment_7_2 )*
            {
             before(grammarAccess.getElkEdgeAccess().getPropertiesAssignment_7_2()); 
            // InternalElkGraph.g:2441:2: ( rule__ElkEdge__PropertiesAssignment_7_2 )*
            loop31:
            do {
                int alt31=2;
                int LA31_0 = input.LA(1);

                if ( (LA31_0==RULE_ID) ) {
                    alt31=1;
                }


                switch (alt31) {
            	case 1 :
            	    // InternalElkGraph.g:2441:3: rule__ElkEdge__PropertiesAssignment_7_2
            	    {
            	    pushFollow(FOLLOW_5);
            	    rule__ElkEdge__PropertiesAssignment_7_2();

            	    state._fsp--;


            	    }
            	    break;

            	default :
            	    break loop31;
                }
            } while (true);

             after(grammarAccess.getElkEdgeAccess().getPropertiesAssignment_7_2()); 

            }


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ElkEdge__Group_7__2__Impl"


    // $ANTLR start "rule__ElkEdge__Group_7__3"
    // InternalElkGraph.g:2449:1: rule__ElkEdge__Group_7__3 : rule__ElkEdge__Group_7__3__Impl rule__ElkEdge__Group_7__4 ;
    public final void rule__ElkEdge__Group_7__3() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:2453:1: ( rule__ElkEdge__Group_7__3__Impl rule__ElkEdge__Group_7__4 )
            // InternalElkGraph.g:2454:2: rule__ElkEdge__Group_7__3__Impl rule__ElkEdge__Group_7__4
            {
            pushFollow(FOLLOW_12);
            rule__ElkEdge__Group_7__3__Impl();

            state._fsp--;

            pushFollow(FOLLOW_2);
            rule__ElkEdge__Group_7__4();

            state._fsp--;


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ElkEdge__Group_7__3"


    // $ANTLR start "rule__ElkEdge__Group_7__3__Impl"
    // InternalElkGraph.g:2461:1: rule__ElkEdge__Group_7__3__Impl : ( ( rule__ElkEdge__LabelsAssignment_7_3 )* ) ;
    public final void rule__ElkEdge__Group_7__3__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:2465:1: ( ( ( rule__ElkEdge__LabelsAssignment_7_3 )* ) )
            // InternalElkGraph.g:2466:1: ( ( rule__ElkEdge__LabelsAssignment_7_3 )* )
            {
            // InternalElkGraph.g:2466:1: ( ( rule__ElkEdge__LabelsAssignment_7_3 )* )
            // InternalElkGraph.g:2467:2: ( rule__ElkEdge__LabelsAssignment_7_3 )*
            {
             before(grammarAccess.getElkEdgeAccess().getLabelsAssignment_7_3()); 
            // InternalElkGraph.g:2468:2: ( rule__ElkEdge__LabelsAssignment_7_3 )*
            loop32:
            do {
                int alt32=2;
                int LA32_0 = input.LA(1);

                if ( (LA32_0==20) ) {
                    alt32=1;
                }


                switch (alt32) {
            	case 1 :
            	    // InternalElkGraph.g:2468:3: rule__ElkEdge__LabelsAssignment_7_3
            	    {
            	    pushFollow(FOLLOW_13);
            	    rule__ElkEdge__LabelsAssignment_7_3();

            	    state._fsp--;


            	    }
            	    break;

            	default :
            	    break loop32;
                }
            } while (true);

             after(grammarAccess.getElkEdgeAccess().getLabelsAssignment_7_3()); 

            }


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ElkEdge__Group_7__3__Impl"


    // $ANTLR start "rule__ElkEdge__Group_7__4"
    // InternalElkGraph.g:2476:1: rule__ElkEdge__Group_7__4 : rule__ElkEdge__Group_7__4__Impl ;
    public final void rule__ElkEdge__Group_7__4() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:2480:1: ( rule__ElkEdge__Group_7__4__Impl )
            // InternalElkGraph.g:2481:2: rule__ElkEdge__Group_7__4__Impl
            {
            pushFollow(FOLLOW_2);
            rule__ElkEdge__Group_7__4__Impl();

            state._fsp--;


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ElkEdge__Group_7__4"


    // $ANTLR start "rule__ElkEdge__Group_7__4__Impl"
    // InternalElkGraph.g:2487:1: rule__ElkEdge__Group_7__4__Impl : ( '}' ) ;
    public final void rule__ElkEdge__Group_7__4__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:2491:1: ( ( '}' ) )
            // InternalElkGraph.g:2492:1: ( '}' )
            {
            // InternalElkGraph.g:2492:1: ( '}' )
            // InternalElkGraph.g:2493:2: '}'
            {
             before(grammarAccess.getElkEdgeAccess().getRightCurlyBracketKeyword_7_4()); 
            match(input,19,FOLLOW_2); 
             after(grammarAccess.getElkEdgeAccess().getRightCurlyBracketKeyword_7_4()); 

            }


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ElkEdge__Group_7__4__Impl"


    // $ANTLR start "rule__EdgeLayout__Group__0"
    // InternalElkGraph.g:2503:1: rule__EdgeLayout__Group__0 : rule__EdgeLayout__Group__0__Impl rule__EdgeLayout__Group__1 ;
    public final void rule__EdgeLayout__Group__0() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:2507:1: ( rule__EdgeLayout__Group__0__Impl rule__EdgeLayout__Group__1 )
            // InternalElkGraph.g:2508:2: rule__EdgeLayout__Group__0__Impl rule__EdgeLayout__Group__1
            {
            pushFollow(FOLLOW_14);
            rule__EdgeLayout__Group__0__Impl();

            state._fsp--;

            pushFollow(FOLLOW_2);
            rule__EdgeLayout__Group__1();

            state._fsp--;


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__EdgeLayout__Group__0"


    // $ANTLR start "rule__EdgeLayout__Group__0__Impl"
    // InternalElkGraph.g:2515:1: rule__EdgeLayout__Group__0__Impl : ( 'layout' ) ;
    public final void rule__EdgeLayout__Group__0__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:2519:1: ( ( 'layout' ) )
            // InternalElkGraph.g:2520:1: ( 'layout' )
            {
            // InternalElkGraph.g:2520:1: ( 'layout' )
            // InternalElkGraph.g:2521:2: 'layout'
            {
             before(grammarAccess.getEdgeLayoutAccess().getLayoutKeyword_0()); 
            match(input,23,FOLLOW_2); 
             after(grammarAccess.getEdgeLayoutAccess().getLayoutKeyword_0()); 

            }


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__EdgeLayout__Group__0__Impl"


    // $ANTLR start "rule__EdgeLayout__Group__1"
    // InternalElkGraph.g:2530:1: rule__EdgeLayout__Group__1 : rule__EdgeLayout__Group__1__Impl rule__EdgeLayout__Group__2 ;
    public final void rule__EdgeLayout__Group__1() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:2534:1: ( rule__EdgeLayout__Group__1__Impl rule__EdgeLayout__Group__2 )
            // InternalElkGraph.g:2535:2: rule__EdgeLayout__Group__1__Impl rule__EdgeLayout__Group__2
            {
            pushFollow(FOLLOW_22);
            rule__EdgeLayout__Group__1__Impl();

            state._fsp--;

            pushFollow(FOLLOW_2);
            rule__EdgeLayout__Group__2();

            state._fsp--;


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__EdgeLayout__Group__1"


    // $ANTLR start "rule__EdgeLayout__Group__1__Impl"
    // InternalElkGraph.g:2542:1: rule__EdgeLayout__Group__1__Impl : ( '[' ) ;
    public final void rule__EdgeLayout__Group__1__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:2546:1: ( ( '[' ) )
            // InternalElkGraph.g:2547:1: ( '[' )
            {
            // InternalElkGraph.g:2547:1: ( '[' )
            // InternalElkGraph.g:2548:2: '['
            {
             before(grammarAccess.getEdgeLayoutAccess().getLeftSquareBracketKeyword_1()); 
            match(input,24,FOLLOW_2); 
             after(grammarAccess.getEdgeLayoutAccess().getLeftSquareBracketKeyword_1()); 

            }


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__EdgeLayout__Group__1__Impl"


    // $ANTLR start "rule__EdgeLayout__Group__2"
    // InternalElkGraph.g:2557:1: rule__EdgeLayout__Group__2 : rule__EdgeLayout__Group__2__Impl rule__EdgeLayout__Group__3 ;
    public final void rule__EdgeLayout__Group__2() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:2561:1: ( rule__EdgeLayout__Group__2__Impl rule__EdgeLayout__Group__3 )
            // InternalElkGraph.g:2562:2: rule__EdgeLayout__Group__2__Impl rule__EdgeLayout__Group__3
            {
            pushFollow(FOLLOW_16);
            rule__EdgeLayout__Group__2__Impl();

            state._fsp--;

            pushFollow(FOLLOW_2);
            rule__EdgeLayout__Group__3();

            state._fsp--;


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__EdgeLayout__Group__2"


    // $ANTLR start "rule__EdgeLayout__Group__2__Impl"
    // InternalElkGraph.g:2569:1: rule__EdgeLayout__Group__2__Impl : ( ( rule__EdgeLayout__Alternatives_2 ) ) ;
    public final void rule__EdgeLayout__Group__2__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:2573:1: ( ( ( rule__EdgeLayout__Alternatives_2 ) ) )
            // InternalElkGraph.g:2574:1: ( ( rule__EdgeLayout__Alternatives_2 ) )
            {
            // InternalElkGraph.g:2574:1: ( ( rule__EdgeLayout__Alternatives_2 ) )
            // InternalElkGraph.g:2575:2: ( rule__EdgeLayout__Alternatives_2 )
            {
             before(grammarAccess.getEdgeLayoutAccess().getAlternatives_2()); 
            // InternalElkGraph.g:2576:2: ( rule__EdgeLayout__Alternatives_2 )
            // InternalElkGraph.g:2576:3: rule__EdgeLayout__Alternatives_2
            {
            pushFollow(FOLLOW_2);
            rule__EdgeLayout__Alternatives_2();

            state._fsp--;


            }

             after(grammarAccess.getEdgeLayoutAccess().getAlternatives_2()); 

            }


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__EdgeLayout__Group__2__Impl"


    // $ANTLR start "rule__EdgeLayout__Group__3"
    // InternalElkGraph.g:2584:1: rule__EdgeLayout__Group__3 : rule__EdgeLayout__Group__3__Impl ;
    public final void rule__EdgeLayout__Group__3() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:2588:1: ( rule__EdgeLayout__Group__3__Impl )
            // InternalElkGraph.g:2589:2: rule__EdgeLayout__Group__3__Impl
            {
            pushFollow(FOLLOW_2);
            rule__EdgeLayout__Group__3__Impl();

            state._fsp--;


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__EdgeLayout__Group__3"


    // $ANTLR start "rule__EdgeLayout__Group__3__Impl"
    // InternalElkGraph.g:2595:1: rule__EdgeLayout__Group__3__Impl : ( ']' ) ;
    public final void rule__EdgeLayout__Group__3__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:2599:1: ( ( ']' ) )
            // InternalElkGraph.g:2600:1: ( ']' )
            {
            // InternalElkGraph.g:2600:1: ( ']' )
            // InternalElkGraph.g:2601:2: ']'
            {
             before(grammarAccess.getEdgeLayoutAccess().getRightSquareBracketKeyword_3()); 
            match(input,25,FOLLOW_2); 
             after(grammarAccess.getEdgeLayoutAccess().getRightSquareBracketKeyword_3()); 

            }


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__EdgeLayout__Group__3__Impl"


    // $ANTLR start "rule__ElkSingleEdgeSection__Group__0"
    // InternalElkGraph.g:2611:1: rule__ElkSingleEdgeSection__Group__0 : rule__ElkSingleEdgeSection__Group__0__Impl rule__ElkSingleEdgeSection__Group__1 ;
    public final void rule__ElkSingleEdgeSection__Group__0() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:2615:1: ( rule__ElkSingleEdgeSection__Group__0__Impl rule__ElkSingleEdgeSection__Group__1 )
            // InternalElkGraph.g:2616:2: rule__ElkSingleEdgeSection__Group__0__Impl rule__ElkSingleEdgeSection__Group__1
            {
            pushFollow(FOLLOW_23);
            rule__ElkSingleEdgeSection__Group__0__Impl();

            state._fsp--;

            pushFollow(FOLLOW_2);
            rule__ElkSingleEdgeSection__Group__1();

            state._fsp--;


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ElkSingleEdgeSection__Group__0"


    // $ANTLR start "rule__ElkSingleEdgeSection__Group__0__Impl"
    // InternalElkGraph.g:2623:1: rule__ElkSingleEdgeSection__Group__0__Impl : ( () ) ;
    public final void rule__ElkSingleEdgeSection__Group__0__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:2627:1: ( ( () ) )
            // InternalElkGraph.g:2628:1: ( () )
            {
            // InternalElkGraph.g:2628:1: ( () )
            // InternalElkGraph.g:2629:2: ()
            {
             before(grammarAccess.getElkSingleEdgeSectionAccess().getElkEdgeSectionAction_0()); 
            // InternalElkGraph.g:2630:2: ()
            // InternalElkGraph.g:2630:3: 
            {
            }

             after(grammarAccess.getElkSingleEdgeSectionAccess().getElkEdgeSectionAction_0()); 

            }


            }

        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ElkSingleEdgeSection__Group__0__Impl"


    // $ANTLR start "rule__ElkSingleEdgeSection__Group__1"
    // InternalElkGraph.g:2638:1: rule__ElkSingleEdgeSection__Group__1 : rule__ElkSingleEdgeSection__Group__1__Impl ;
    public final void rule__ElkSingleEdgeSection__Group__1() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:2642:1: ( rule__ElkSingleEdgeSection__Group__1__Impl )
            // InternalElkGraph.g:2643:2: rule__ElkSingleEdgeSection__Group__1__Impl
            {
            pushFollow(FOLLOW_2);
            rule__ElkSingleEdgeSection__Group__1__Impl();

            state._fsp--;


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ElkSingleEdgeSection__Group__1"


    // $ANTLR start "rule__ElkSingleEdgeSection__Group__1__Impl"
    // InternalElkGraph.g:2649:1: rule__ElkSingleEdgeSection__Group__1__Impl : ( ( rule__ElkSingleEdgeSection__Group_1__0 ) ) ;
    public final void rule__ElkSingleEdgeSection__Group__1__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:2653:1: ( ( ( rule__ElkSingleEdgeSection__Group_1__0 ) ) )
            // InternalElkGraph.g:2654:1: ( ( rule__ElkSingleEdgeSection__Group_1__0 ) )
            {
            // InternalElkGraph.g:2654:1: ( ( rule__ElkSingleEdgeSection__Group_1__0 ) )
            // InternalElkGraph.g:2655:2: ( rule__ElkSingleEdgeSection__Group_1__0 )
            {
             before(grammarAccess.getElkSingleEdgeSectionAccess().getGroup_1()); 
            // InternalElkGraph.g:2656:2: ( rule__ElkSingleEdgeSection__Group_1__0 )
            // InternalElkGraph.g:2656:3: rule__ElkSingleEdgeSection__Group_1__0
            {
            pushFollow(FOLLOW_2);
            rule__ElkSingleEdgeSection__Group_1__0();

            state._fsp--;


            }

             after(grammarAccess.getElkSingleEdgeSectionAccess().getGroup_1()); 

            }


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ElkSingleEdgeSection__Group__1__Impl"


    // $ANTLR start "rule__ElkSingleEdgeSection__Group_1__0"
    // InternalElkGraph.g:2665:1: rule__ElkSingleEdgeSection__Group_1__0 : rule__ElkSingleEdgeSection__Group_1__0__Impl rule__ElkSingleEdgeSection__Group_1__1 ;
    public final void rule__ElkSingleEdgeSection__Group_1__0() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:2669:1: ( rule__ElkSingleEdgeSection__Group_1__0__Impl rule__ElkSingleEdgeSection__Group_1__1 )
            // InternalElkGraph.g:2670:2: rule__ElkSingleEdgeSection__Group_1__0__Impl rule__ElkSingleEdgeSection__Group_1__1
            {
            pushFollow(FOLLOW_24);
            rule__ElkSingleEdgeSection__Group_1__0__Impl();

            state._fsp--;

            pushFollow(FOLLOW_2);
            rule__ElkSingleEdgeSection__Group_1__1();

            state._fsp--;


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ElkSingleEdgeSection__Group_1__0"


    // $ANTLR start "rule__ElkSingleEdgeSection__Group_1__0__Impl"
    // InternalElkGraph.g:2677:1: rule__ElkSingleEdgeSection__Group_1__0__Impl : ( ( rule__ElkSingleEdgeSection__UnorderedGroup_1_0 ) ) ;
    public final void rule__ElkSingleEdgeSection__Group_1__0__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:2681:1: ( ( ( rule__ElkSingleEdgeSection__UnorderedGroup_1_0 ) ) )
            // InternalElkGraph.g:2682:1: ( ( rule__ElkSingleEdgeSection__UnorderedGroup_1_0 ) )
            {
            // InternalElkGraph.g:2682:1: ( ( rule__ElkSingleEdgeSection__UnorderedGroup_1_0 ) )
            // InternalElkGraph.g:2683:2: ( rule__ElkSingleEdgeSection__UnorderedGroup_1_0 )
            {
             before(grammarAccess.getElkSingleEdgeSectionAccess().getUnorderedGroup_1_0()); 
            // InternalElkGraph.g:2684:2: ( rule__ElkSingleEdgeSection__UnorderedGroup_1_0 )
            // InternalElkGraph.g:2684:3: rule__ElkSingleEdgeSection__UnorderedGroup_1_0
            {
            pushFollow(FOLLOW_2);
            rule__ElkSingleEdgeSection__UnorderedGroup_1_0();

            state._fsp--;


            }

             after(grammarAccess.getElkSingleEdgeSectionAccess().getUnorderedGroup_1_0()); 

            }


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ElkSingleEdgeSection__Group_1__0__Impl"


    // $ANTLR start "rule__ElkSingleEdgeSection__Group_1__1"
    // InternalElkGraph.g:2692:1: rule__ElkSingleEdgeSection__Group_1__1 : rule__ElkSingleEdgeSection__Group_1__1__Impl rule__ElkSingleEdgeSection__Group_1__2 ;
    public final void rule__ElkSingleEdgeSection__Group_1__1() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:2696:1: ( rule__ElkSingleEdgeSection__Group_1__1__Impl rule__ElkSingleEdgeSection__Group_1__2 )
            // InternalElkGraph.g:2697:2: rule__ElkSingleEdgeSection__Group_1__1__Impl rule__ElkSingleEdgeSection__Group_1__2
            {
            pushFollow(FOLLOW_24);
            rule__ElkSingleEdgeSection__Group_1__1__Impl();

            state._fsp--;

            pushFollow(FOLLOW_2);
            rule__ElkSingleEdgeSection__Group_1__2();

            state._fsp--;


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ElkSingleEdgeSection__Group_1__1"


    // $ANTLR start "rule__ElkSingleEdgeSection__Group_1__1__Impl"
    // InternalElkGraph.g:2704:1: rule__ElkSingleEdgeSection__Group_1__1__Impl : ( ( rule__ElkSingleEdgeSection__Group_1_1__0 )? ) ;
    public final void rule__ElkSingleEdgeSection__Group_1__1__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:2708:1: ( ( ( rule__ElkSingleEdgeSection__Group_1_1__0 )? ) )
            // InternalElkGraph.g:2709:1: ( ( rule__ElkSingleEdgeSection__Group_1_1__0 )? )
            {
            // InternalElkGraph.g:2709:1: ( ( rule__ElkSingleEdgeSection__Group_1_1__0 )? )
            // InternalElkGraph.g:2710:2: ( rule__ElkSingleEdgeSection__Group_1_1__0 )?
            {
             before(grammarAccess.getElkSingleEdgeSectionAccess().getGroup_1_1()); 
            // InternalElkGraph.g:2711:2: ( rule__ElkSingleEdgeSection__Group_1_1__0 )?
            int alt33=2;
            int LA33_0 = input.LA(1);

            if ( (LA33_0==35) ) {
                alt33=1;
            }
            switch (alt33) {
                case 1 :
                    // InternalElkGraph.g:2711:3: rule__ElkSingleEdgeSection__Group_1_1__0
                    {
                    pushFollow(FOLLOW_2);
                    rule__ElkSingleEdgeSection__Group_1_1__0();

                    state._fsp--;


                    }
                    break;

            }

             after(grammarAccess.getElkSingleEdgeSectionAccess().getGroup_1_1()); 

            }


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ElkSingleEdgeSection__Group_1__1__Impl"


    // $ANTLR start "rule__ElkSingleEdgeSection__Group_1__2"
    // InternalElkGraph.g:2719:1: rule__ElkSingleEdgeSection__Group_1__2 : rule__ElkSingleEdgeSection__Group_1__2__Impl ;
    public final void rule__ElkSingleEdgeSection__Group_1__2() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:2723:1: ( rule__ElkSingleEdgeSection__Group_1__2__Impl )
            // InternalElkGraph.g:2724:2: rule__ElkSingleEdgeSection__Group_1__2__Impl
            {
            pushFollow(FOLLOW_2);
            rule__ElkSingleEdgeSection__Group_1__2__Impl();

            state._fsp--;


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ElkSingleEdgeSection__Group_1__2"


    // $ANTLR start "rule__ElkSingleEdgeSection__Group_1__2__Impl"
    // InternalElkGraph.g:2730:1: rule__ElkSingleEdgeSection__Group_1__2__Impl : ( ( rule__ElkSingleEdgeSection__PropertiesAssignment_1_2 )* ) ;
    public final void rule__ElkSingleEdgeSection__Group_1__2__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:2734:1: ( ( ( rule__ElkSingleEdgeSection__PropertiesAssignment_1_2 )* ) )
            // InternalElkGraph.g:2735:1: ( ( rule__ElkSingleEdgeSection__PropertiesAssignment_1_2 )* )
            {
            // InternalElkGraph.g:2735:1: ( ( rule__ElkSingleEdgeSection__PropertiesAssignment_1_2 )* )
            // InternalElkGraph.g:2736:2: ( rule__ElkSingleEdgeSection__PropertiesAssignment_1_2 )*
            {
             before(grammarAccess.getElkSingleEdgeSectionAccess().getPropertiesAssignment_1_2()); 
            // InternalElkGraph.g:2737:2: ( rule__ElkSingleEdgeSection__PropertiesAssignment_1_2 )*
            loop34:
            do {
                int alt34=2;
                int LA34_0 = input.LA(1);

                if ( (LA34_0==RULE_ID) ) {
                    alt34=1;
                }


                switch (alt34) {
            	case 1 :
            	    // InternalElkGraph.g:2737:3: rule__ElkSingleEdgeSection__PropertiesAssignment_1_2
            	    {
            	    pushFollow(FOLLOW_5);
            	    rule__ElkSingleEdgeSection__PropertiesAssignment_1_2();

            	    state._fsp--;


            	    }
            	    break;

            	default :
            	    break loop34;
                }
            } while (true);

             after(grammarAccess.getElkSingleEdgeSectionAccess().getPropertiesAssignment_1_2()); 

            }


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ElkSingleEdgeSection__Group_1__2__Impl"


    // $ANTLR start "rule__ElkSingleEdgeSection__Group_1_0_0__0"
    // InternalElkGraph.g:2746:1: rule__ElkSingleEdgeSection__Group_1_0_0__0 : rule__ElkSingleEdgeSection__Group_1_0_0__0__Impl rule__ElkSingleEdgeSection__Group_1_0_0__1 ;
    public final void rule__ElkSingleEdgeSection__Group_1_0_0__0() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:2750:1: ( rule__ElkSingleEdgeSection__Group_1_0_0__0__Impl rule__ElkSingleEdgeSection__Group_1_0_0__1 )
            // InternalElkGraph.g:2751:2: rule__ElkSingleEdgeSection__Group_1_0_0__0__Impl rule__ElkSingleEdgeSection__Group_1_0_0__1
            {
            pushFollow(FOLLOW_11);
            rule__ElkSingleEdgeSection__Group_1_0_0__0__Impl();

            state._fsp--;

            pushFollow(FOLLOW_2);
            rule__ElkSingleEdgeSection__Group_1_0_0__1();

            state._fsp--;


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ElkSingleEdgeSection__Group_1_0_0__0"


    // $ANTLR start "rule__ElkSingleEdgeSection__Group_1_0_0__0__Impl"
    // InternalElkGraph.g:2758:1: rule__ElkSingleEdgeSection__Group_1_0_0__0__Impl : ( 'incoming' ) ;
    public final void rule__ElkSingleEdgeSection__Group_1_0_0__0__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:2762:1: ( ( 'incoming' ) )
            // InternalElkGraph.g:2763:1: ( 'incoming' )
            {
            // InternalElkGraph.g:2763:1: ( 'incoming' )
            // InternalElkGraph.g:2764:2: 'incoming'
            {
             before(grammarAccess.getElkSingleEdgeSectionAccess().getIncomingKeyword_1_0_0_0()); 
            match(input,31,FOLLOW_2); 
             after(grammarAccess.getElkSingleEdgeSectionAccess().getIncomingKeyword_1_0_0_0()); 

            }


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ElkSingleEdgeSection__Group_1_0_0__0__Impl"


    // $ANTLR start "rule__ElkSingleEdgeSection__Group_1_0_0__1"
    // InternalElkGraph.g:2773:1: rule__ElkSingleEdgeSection__Group_1_0_0__1 : rule__ElkSingleEdgeSection__Group_1_0_0__1__Impl rule__ElkSingleEdgeSection__Group_1_0_0__2 ;
    public final void rule__ElkSingleEdgeSection__Group_1_0_0__1() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:2777:1: ( rule__ElkSingleEdgeSection__Group_1_0_0__1__Impl rule__ElkSingleEdgeSection__Group_1_0_0__2 )
            // InternalElkGraph.g:2778:2: rule__ElkSingleEdgeSection__Group_1_0_0__1__Impl rule__ElkSingleEdgeSection__Group_1_0_0__2
            {
            pushFollow(FOLLOW_7);
            rule__ElkSingleEdgeSection__Group_1_0_0__1__Impl();

            state._fsp--;

            pushFollow(FOLLOW_2);
            rule__ElkSingleEdgeSection__Group_1_0_0__2();

            state._fsp--;


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ElkSingleEdgeSection__Group_1_0_0__1"


    // $ANTLR start "rule__ElkSingleEdgeSection__Group_1_0_0__1__Impl"
    // InternalElkGraph.g:2785:1: rule__ElkSingleEdgeSection__Group_1_0_0__1__Impl : ( ':' ) ;
    public final void rule__ElkSingleEdgeSection__Group_1_0_0__1__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:2789:1: ( ( ':' ) )
            // InternalElkGraph.g:2790:1: ( ':' )
            {
            // InternalElkGraph.g:2790:1: ( ':' )
            // InternalElkGraph.g:2791:2: ':'
            {
             before(grammarAccess.getElkSingleEdgeSectionAccess().getColonKeyword_1_0_0_1()); 
            match(input,21,FOLLOW_2); 
             after(grammarAccess.getElkSingleEdgeSectionAccess().getColonKeyword_1_0_0_1()); 

            }


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ElkSingleEdgeSection__Group_1_0_0__1__Impl"


    // $ANTLR start "rule__ElkSingleEdgeSection__Group_1_0_0__2"
    // InternalElkGraph.g:2800:1: rule__ElkSingleEdgeSection__Group_1_0_0__2 : rule__ElkSingleEdgeSection__Group_1_0_0__2__Impl ;
    public final void rule__ElkSingleEdgeSection__Group_1_0_0__2() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:2804:1: ( rule__ElkSingleEdgeSection__Group_1_0_0__2__Impl )
            // InternalElkGraph.g:2805:2: rule__ElkSingleEdgeSection__Group_1_0_0__2__Impl
            {
            pushFollow(FOLLOW_2);
            rule__ElkSingleEdgeSection__Group_1_0_0__2__Impl();

            state._fsp--;


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ElkSingleEdgeSection__Group_1_0_0__2"


    // $ANTLR start "rule__ElkSingleEdgeSection__Group_1_0_0__2__Impl"
    // InternalElkGraph.g:2811:1: rule__ElkSingleEdgeSection__Group_1_0_0__2__Impl : ( ( rule__ElkSingleEdgeSection__IncomingShapeAssignment_1_0_0_2 ) ) ;
    public final void rule__ElkSingleEdgeSection__Group_1_0_0__2__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:2815:1: ( ( ( rule__ElkSingleEdgeSection__IncomingShapeAssignment_1_0_0_2 ) ) )
            // InternalElkGraph.g:2816:1: ( ( rule__ElkSingleEdgeSection__IncomingShapeAssignment_1_0_0_2 ) )
            {
            // InternalElkGraph.g:2816:1: ( ( rule__ElkSingleEdgeSection__IncomingShapeAssignment_1_0_0_2 ) )
            // InternalElkGraph.g:2817:2: ( rule__ElkSingleEdgeSection__IncomingShapeAssignment_1_0_0_2 )
            {
             before(grammarAccess.getElkSingleEdgeSectionAccess().getIncomingShapeAssignment_1_0_0_2()); 
            // InternalElkGraph.g:2818:2: ( rule__ElkSingleEdgeSection__IncomingShapeAssignment_1_0_0_2 )
            // InternalElkGraph.g:2818:3: rule__ElkSingleEdgeSection__IncomingShapeAssignment_1_0_0_2
            {
            pushFollow(FOLLOW_2);
            rule__ElkSingleEdgeSection__IncomingShapeAssignment_1_0_0_2();

            state._fsp--;


            }

             after(grammarAccess.getElkSingleEdgeSectionAccess().getIncomingShapeAssignment_1_0_0_2()); 

            }


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ElkSingleEdgeSection__Group_1_0_0__2__Impl"


    // $ANTLR start "rule__ElkSingleEdgeSection__Group_1_0_1__0"
    // InternalElkGraph.g:2827:1: rule__ElkSingleEdgeSection__Group_1_0_1__0 : rule__ElkSingleEdgeSection__Group_1_0_1__0__Impl rule__ElkSingleEdgeSection__Group_1_0_1__1 ;
    public final void rule__ElkSingleEdgeSection__Group_1_0_1__0() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:2831:1: ( rule__ElkSingleEdgeSection__Group_1_0_1__0__Impl rule__ElkSingleEdgeSection__Group_1_0_1__1 )
            // InternalElkGraph.g:2832:2: rule__ElkSingleEdgeSection__Group_1_0_1__0__Impl rule__ElkSingleEdgeSection__Group_1_0_1__1
            {
            pushFollow(FOLLOW_11);
            rule__ElkSingleEdgeSection__Group_1_0_1__0__Impl();

            state._fsp--;

            pushFollow(FOLLOW_2);
            rule__ElkSingleEdgeSection__Group_1_0_1__1();

            state._fsp--;


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ElkSingleEdgeSection__Group_1_0_1__0"


    // $ANTLR start "rule__ElkSingleEdgeSection__Group_1_0_1__0__Impl"
    // InternalElkGraph.g:2839:1: rule__ElkSingleEdgeSection__Group_1_0_1__0__Impl : ( 'outgoing' ) ;
    public final void rule__ElkSingleEdgeSection__Group_1_0_1__0__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:2843:1: ( ( 'outgoing' ) )
            // InternalElkGraph.g:2844:1: ( 'outgoing' )
            {
            // InternalElkGraph.g:2844:1: ( 'outgoing' )
            // InternalElkGraph.g:2845:2: 'outgoing'
            {
             before(grammarAccess.getElkSingleEdgeSectionAccess().getOutgoingKeyword_1_0_1_0()); 
            match(input,32,FOLLOW_2); 
             after(grammarAccess.getElkSingleEdgeSectionAccess().getOutgoingKeyword_1_0_1_0()); 

            }


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ElkSingleEdgeSection__Group_1_0_1__0__Impl"


    // $ANTLR start "rule__ElkSingleEdgeSection__Group_1_0_1__1"
    // InternalElkGraph.g:2854:1: rule__ElkSingleEdgeSection__Group_1_0_1__1 : rule__ElkSingleEdgeSection__Group_1_0_1__1__Impl rule__ElkSingleEdgeSection__Group_1_0_1__2 ;
    public final void rule__ElkSingleEdgeSection__Group_1_0_1__1() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:2858:1: ( rule__ElkSingleEdgeSection__Group_1_0_1__1__Impl rule__ElkSingleEdgeSection__Group_1_0_1__2 )
            // InternalElkGraph.g:2859:2: rule__ElkSingleEdgeSection__Group_1_0_1__1__Impl rule__ElkSingleEdgeSection__Group_1_0_1__2
            {
            pushFollow(FOLLOW_7);
            rule__ElkSingleEdgeSection__Group_1_0_1__1__Impl();

            state._fsp--;

            pushFollow(FOLLOW_2);
            rule__ElkSingleEdgeSection__Group_1_0_1__2();

            state._fsp--;


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ElkSingleEdgeSection__Group_1_0_1__1"


    // $ANTLR start "rule__ElkSingleEdgeSection__Group_1_0_1__1__Impl"
    // InternalElkGraph.g:2866:1: rule__ElkSingleEdgeSection__Group_1_0_1__1__Impl : ( ':' ) ;
    public final void rule__ElkSingleEdgeSection__Group_1_0_1__1__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:2870:1: ( ( ':' ) )
            // InternalElkGraph.g:2871:1: ( ':' )
            {
            // InternalElkGraph.g:2871:1: ( ':' )
            // InternalElkGraph.g:2872:2: ':'
            {
             before(grammarAccess.getElkSingleEdgeSectionAccess().getColonKeyword_1_0_1_1()); 
            match(input,21,FOLLOW_2); 
             after(grammarAccess.getElkSingleEdgeSectionAccess().getColonKeyword_1_0_1_1()); 

            }


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ElkSingleEdgeSection__Group_1_0_1__1__Impl"


    // $ANTLR start "rule__ElkSingleEdgeSection__Group_1_0_1__2"
    // InternalElkGraph.g:2881:1: rule__ElkSingleEdgeSection__Group_1_0_1__2 : rule__ElkSingleEdgeSection__Group_1_0_1__2__Impl ;
    public final void rule__ElkSingleEdgeSection__Group_1_0_1__2() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:2885:1: ( rule__ElkSingleEdgeSection__Group_1_0_1__2__Impl )
            // InternalElkGraph.g:2886:2: rule__ElkSingleEdgeSection__Group_1_0_1__2__Impl
            {
            pushFollow(FOLLOW_2);
            rule__ElkSingleEdgeSection__Group_1_0_1__2__Impl();

            state._fsp--;


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ElkSingleEdgeSection__Group_1_0_1__2"


    // $ANTLR start "rule__ElkSingleEdgeSection__Group_1_0_1__2__Impl"
    // InternalElkGraph.g:2892:1: rule__ElkSingleEdgeSection__Group_1_0_1__2__Impl : ( ( rule__ElkSingleEdgeSection__OutgoingShapeAssignment_1_0_1_2 ) ) ;
    public final void rule__ElkSingleEdgeSection__Group_1_0_1__2__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:2896:1: ( ( ( rule__ElkSingleEdgeSection__OutgoingShapeAssignment_1_0_1_2 ) ) )
            // InternalElkGraph.g:2897:1: ( ( rule__ElkSingleEdgeSection__OutgoingShapeAssignment_1_0_1_2 ) )
            {
            // InternalElkGraph.g:2897:1: ( ( rule__ElkSingleEdgeSection__OutgoingShapeAssignment_1_0_1_2 ) )
            // InternalElkGraph.g:2898:2: ( rule__ElkSingleEdgeSection__OutgoingShapeAssignment_1_0_1_2 )
            {
             before(grammarAccess.getElkSingleEdgeSectionAccess().getOutgoingShapeAssignment_1_0_1_2()); 
            // InternalElkGraph.g:2899:2: ( rule__ElkSingleEdgeSection__OutgoingShapeAssignment_1_0_1_2 )
            // InternalElkGraph.g:2899:3: rule__ElkSingleEdgeSection__OutgoingShapeAssignment_1_0_1_2
            {
            pushFollow(FOLLOW_2);
            rule__ElkSingleEdgeSection__OutgoingShapeAssignment_1_0_1_2();

            state._fsp--;


            }

             after(grammarAccess.getElkSingleEdgeSectionAccess().getOutgoingShapeAssignment_1_0_1_2()); 

            }


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ElkSingleEdgeSection__Group_1_0_1__2__Impl"


    // $ANTLR start "rule__ElkSingleEdgeSection__Group_1_0_2__0"
    // InternalElkGraph.g:2908:1: rule__ElkSingleEdgeSection__Group_1_0_2__0 : rule__ElkSingleEdgeSection__Group_1_0_2__0__Impl rule__ElkSingleEdgeSection__Group_1_0_2__1 ;
    public final void rule__ElkSingleEdgeSection__Group_1_0_2__0() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:2912:1: ( rule__ElkSingleEdgeSection__Group_1_0_2__0__Impl rule__ElkSingleEdgeSection__Group_1_0_2__1 )
            // InternalElkGraph.g:2913:2: rule__ElkSingleEdgeSection__Group_1_0_2__0__Impl rule__ElkSingleEdgeSection__Group_1_0_2__1
            {
            pushFollow(FOLLOW_11);
            rule__ElkSingleEdgeSection__Group_1_0_2__0__Impl();

            state._fsp--;

            pushFollow(FOLLOW_2);
            rule__ElkSingleEdgeSection__Group_1_0_2__1();

            state._fsp--;


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ElkSingleEdgeSection__Group_1_0_2__0"


    // $ANTLR start "rule__ElkSingleEdgeSection__Group_1_0_2__0__Impl"
    // InternalElkGraph.g:2920:1: rule__ElkSingleEdgeSection__Group_1_0_2__0__Impl : ( 'start' ) ;
    public final void rule__ElkSingleEdgeSection__Group_1_0_2__0__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:2924:1: ( ( 'start' ) )
            // InternalElkGraph.g:2925:1: ( 'start' )
            {
            // InternalElkGraph.g:2925:1: ( 'start' )
            // InternalElkGraph.g:2926:2: 'start'
            {
             before(grammarAccess.getElkSingleEdgeSectionAccess().getStartKeyword_1_0_2_0()); 
            match(input,33,FOLLOW_2); 
             after(grammarAccess.getElkSingleEdgeSectionAccess().getStartKeyword_1_0_2_0()); 

            }


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ElkSingleEdgeSection__Group_1_0_2__0__Impl"


    // $ANTLR start "rule__ElkSingleEdgeSection__Group_1_0_2__1"
    // InternalElkGraph.g:2935:1: rule__ElkSingleEdgeSection__Group_1_0_2__1 : rule__ElkSingleEdgeSection__Group_1_0_2__1__Impl rule__ElkSingleEdgeSection__Group_1_0_2__2 ;
    public final void rule__ElkSingleEdgeSection__Group_1_0_2__1() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:2939:1: ( rule__ElkSingleEdgeSection__Group_1_0_2__1__Impl rule__ElkSingleEdgeSection__Group_1_0_2__2 )
            // InternalElkGraph.g:2940:2: rule__ElkSingleEdgeSection__Group_1_0_2__1__Impl rule__ElkSingleEdgeSection__Group_1_0_2__2
            {
            pushFollow(FOLLOW_17);
            rule__ElkSingleEdgeSection__Group_1_0_2__1__Impl();

            state._fsp--;

            pushFollow(FOLLOW_2);
            rule__ElkSingleEdgeSection__Group_1_0_2__2();

            state._fsp--;


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ElkSingleEdgeSection__Group_1_0_2__1"


    // $ANTLR start "rule__ElkSingleEdgeSection__Group_1_0_2__1__Impl"
    // InternalElkGraph.g:2947:1: rule__ElkSingleEdgeSection__Group_1_0_2__1__Impl : ( ':' ) ;
    public final void rule__ElkSingleEdgeSection__Group_1_0_2__1__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:2951:1: ( ( ':' ) )
            // InternalElkGraph.g:2952:1: ( ':' )
            {
            // InternalElkGraph.g:2952:1: ( ':' )
            // InternalElkGraph.g:2953:2: ':'
            {
             before(grammarAccess.getElkSingleEdgeSectionAccess().getColonKeyword_1_0_2_1()); 
            match(input,21,FOLLOW_2); 
             after(grammarAccess.getElkSingleEdgeSectionAccess().getColonKeyword_1_0_2_1()); 

            }


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ElkSingleEdgeSection__Group_1_0_2__1__Impl"


    // $ANTLR start "rule__ElkSingleEdgeSection__Group_1_0_2__2"
    // InternalElkGraph.g:2962:1: rule__ElkSingleEdgeSection__Group_1_0_2__2 : rule__ElkSingleEdgeSection__Group_1_0_2__2__Impl rule__ElkSingleEdgeSection__Group_1_0_2__3 ;
    public final void rule__ElkSingleEdgeSection__Group_1_0_2__2() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:2966:1: ( rule__ElkSingleEdgeSection__Group_1_0_2__2__Impl rule__ElkSingleEdgeSection__Group_1_0_2__3 )
            // InternalElkGraph.g:2967:2: rule__ElkSingleEdgeSection__Group_1_0_2__2__Impl rule__ElkSingleEdgeSection__Group_1_0_2__3
            {
            pushFollow(FOLLOW_18);
            rule__ElkSingleEdgeSection__Group_1_0_2__2__Impl();

            state._fsp--;

            pushFollow(FOLLOW_2);
            rule__ElkSingleEdgeSection__Group_1_0_2__3();

            state._fsp--;


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ElkSingleEdgeSection__Group_1_0_2__2"


    // $ANTLR start "rule__ElkSingleEdgeSection__Group_1_0_2__2__Impl"
    // InternalElkGraph.g:2974:1: rule__ElkSingleEdgeSection__Group_1_0_2__2__Impl : ( ( rule__ElkSingleEdgeSection__StartXAssignment_1_0_2_2 ) ) ;
    public final void rule__ElkSingleEdgeSection__Group_1_0_2__2__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:2978:1: ( ( ( rule__ElkSingleEdgeSection__StartXAssignment_1_0_2_2 ) ) )
            // InternalElkGraph.g:2979:1: ( ( rule__ElkSingleEdgeSection__StartXAssignment_1_0_2_2 ) )
            {
            // InternalElkGraph.g:2979:1: ( ( rule__ElkSingleEdgeSection__StartXAssignment_1_0_2_2 ) )
            // InternalElkGraph.g:2980:2: ( rule__ElkSingleEdgeSection__StartXAssignment_1_0_2_2 )
            {
             before(grammarAccess.getElkSingleEdgeSectionAccess().getStartXAssignment_1_0_2_2()); 
            // InternalElkGraph.g:2981:2: ( rule__ElkSingleEdgeSection__StartXAssignment_1_0_2_2 )
            // InternalElkGraph.g:2981:3: rule__ElkSingleEdgeSection__StartXAssignment_1_0_2_2
            {
            pushFollow(FOLLOW_2);
            rule__ElkSingleEdgeSection__StartXAssignment_1_0_2_2();

            state._fsp--;


            }

             after(grammarAccess.getElkSingleEdgeSectionAccess().getStartXAssignment_1_0_2_2()); 

            }


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ElkSingleEdgeSection__Group_1_0_2__2__Impl"


    // $ANTLR start "rule__ElkSingleEdgeSection__Group_1_0_2__3"
    // InternalElkGraph.g:2989:1: rule__ElkSingleEdgeSection__Group_1_0_2__3 : rule__ElkSingleEdgeSection__Group_1_0_2__3__Impl rule__ElkSingleEdgeSection__Group_1_0_2__4 ;
    public final void rule__ElkSingleEdgeSection__Group_1_0_2__3() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:2993:1: ( rule__ElkSingleEdgeSection__Group_1_0_2__3__Impl rule__ElkSingleEdgeSection__Group_1_0_2__4 )
            // InternalElkGraph.g:2994:2: rule__ElkSingleEdgeSection__Group_1_0_2__3__Impl rule__ElkSingleEdgeSection__Group_1_0_2__4
            {
            pushFollow(FOLLOW_17);
            rule__ElkSingleEdgeSection__Group_1_0_2__3__Impl();

            state._fsp--;

            pushFollow(FOLLOW_2);
            rule__ElkSingleEdgeSection__Group_1_0_2__4();

            state._fsp--;


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ElkSingleEdgeSection__Group_1_0_2__3"


    // $ANTLR start "rule__ElkSingleEdgeSection__Group_1_0_2__3__Impl"
    // InternalElkGraph.g:3001:1: rule__ElkSingleEdgeSection__Group_1_0_2__3__Impl : ( ',' ) ;
    public final void rule__ElkSingleEdgeSection__Group_1_0_2__3__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:3005:1: ( ( ',' ) )
            // InternalElkGraph.g:3006:1: ( ',' )
            {
            // InternalElkGraph.g:3006:1: ( ',' )
            // InternalElkGraph.g:3007:2: ','
            {
             before(grammarAccess.getElkSingleEdgeSectionAccess().getCommaKeyword_1_0_2_3()); 
            match(input,27,FOLLOW_2); 
             after(grammarAccess.getElkSingleEdgeSectionAccess().getCommaKeyword_1_0_2_3()); 

            }


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ElkSingleEdgeSection__Group_1_0_2__3__Impl"


    // $ANTLR start "rule__ElkSingleEdgeSection__Group_1_0_2__4"
    // InternalElkGraph.g:3016:1: rule__ElkSingleEdgeSection__Group_1_0_2__4 : rule__ElkSingleEdgeSection__Group_1_0_2__4__Impl ;
    public final void rule__ElkSingleEdgeSection__Group_1_0_2__4() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:3020:1: ( rule__ElkSingleEdgeSection__Group_1_0_2__4__Impl )
            // InternalElkGraph.g:3021:2: rule__ElkSingleEdgeSection__Group_1_0_2__4__Impl
            {
            pushFollow(FOLLOW_2);
            rule__ElkSingleEdgeSection__Group_1_0_2__4__Impl();

            state._fsp--;


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ElkSingleEdgeSection__Group_1_0_2__4"


    // $ANTLR start "rule__ElkSingleEdgeSection__Group_1_0_2__4__Impl"
    // InternalElkGraph.g:3027:1: rule__ElkSingleEdgeSection__Group_1_0_2__4__Impl : ( ( rule__ElkSingleEdgeSection__StartYAssignment_1_0_2_4 ) ) ;
    public final void rule__ElkSingleEdgeSection__Group_1_0_2__4__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:3031:1: ( ( ( rule__ElkSingleEdgeSection__StartYAssignment_1_0_2_4 ) ) )
            // InternalElkGraph.g:3032:1: ( ( rule__ElkSingleEdgeSection__StartYAssignment_1_0_2_4 ) )
            {
            // InternalElkGraph.g:3032:1: ( ( rule__ElkSingleEdgeSection__StartYAssignment_1_0_2_4 ) )
            // InternalElkGraph.g:3033:2: ( rule__ElkSingleEdgeSection__StartYAssignment_1_0_2_4 )
            {
             before(grammarAccess.getElkSingleEdgeSectionAccess().getStartYAssignment_1_0_2_4()); 
            // InternalElkGraph.g:3034:2: ( rule__ElkSingleEdgeSection__StartYAssignment_1_0_2_4 )
            // InternalElkGraph.g:3034:3: rule__ElkSingleEdgeSection__StartYAssignment_1_0_2_4
            {
            pushFollow(FOLLOW_2);
            rule__ElkSingleEdgeSection__StartYAssignment_1_0_2_4();

            state._fsp--;


            }

             after(grammarAccess.getElkSingleEdgeSectionAccess().getStartYAssignment_1_0_2_4()); 

            }


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ElkSingleEdgeSection__Group_1_0_2__4__Impl"


    // $ANTLR start "rule__ElkSingleEdgeSection__Group_1_0_3__0"
    // InternalElkGraph.g:3043:1: rule__ElkSingleEdgeSection__Group_1_0_3__0 : rule__ElkSingleEdgeSection__Group_1_0_3__0__Impl rule__ElkSingleEdgeSection__Group_1_0_3__1 ;
    public final void rule__ElkSingleEdgeSection__Group_1_0_3__0() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:3047:1: ( rule__ElkSingleEdgeSection__Group_1_0_3__0__Impl rule__ElkSingleEdgeSection__Group_1_0_3__1 )
            // InternalElkGraph.g:3048:2: rule__ElkSingleEdgeSection__Group_1_0_3__0__Impl rule__ElkSingleEdgeSection__Group_1_0_3__1
            {
            pushFollow(FOLLOW_11);
            rule__ElkSingleEdgeSection__Group_1_0_3__0__Impl();

            state._fsp--;

            pushFollow(FOLLOW_2);
            rule__ElkSingleEdgeSection__Group_1_0_3__1();

            state._fsp--;


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ElkSingleEdgeSection__Group_1_0_3__0"


    // $ANTLR start "rule__ElkSingleEdgeSection__Group_1_0_3__0__Impl"
    // InternalElkGraph.g:3055:1: rule__ElkSingleEdgeSection__Group_1_0_3__0__Impl : ( 'end' ) ;
    public final void rule__ElkSingleEdgeSection__Group_1_0_3__0__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:3059:1: ( ( 'end' ) )
            // InternalElkGraph.g:3060:1: ( 'end' )
            {
            // InternalElkGraph.g:3060:1: ( 'end' )
            // InternalElkGraph.g:3061:2: 'end'
            {
             before(grammarAccess.getElkSingleEdgeSectionAccess().getEndKeyword_1_0_3_0()); 
            match(input,34,FOLLOW_2); 
             after(grammarAccess.getElkSingleEdgeSectionAccess().getEndKeyword_1_0_3_0()); 

            }


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ElkSingleEdgeSection__Group_1_0_3__0__Impl"


    // $ANTLR start "rule__ElkSingleEdgeSection__Group_1_0_3__1"
    // InternalElkGraph.g:3070:1: rule__ElkSingleEdgeSection__Group_1_0_3__1 : rule__ElkSingleEdgeSection__Group_1_0_3__1__Impl rule__ElkSingleEdgeSection__Group_1_0_3__2 ;
    public final void rule__ElkSingleEdgeSection__Group_1_0_3__1() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:3074:1: ( rule__ElkSingleEdgeSection__Group_1_0_3__1__Impl rule__ElkSingleEdgeSection__Group_1_0_3__2 )
            // InternalElkGraph.g:3075:2: rule__ElkSingleEdgeSection__Group_1_0_3__1__Impl rule__ElkSingleEdgeSection__Group_1_0_3__2
            {
            pushFollow(FOLLOW_17);
            rule__ElkSingleEdgeSection__Group_1_0_3__1__Impl();

            state._fsp--;

            pushFollow(FOLLOW_2);
            rule__ElkSingleEdgeSection__Group_1_0_3__2();

            state._fsp--;


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ElkSingleEdgeSection__Group_1_0_3__1"


    // $ANTLR start "rule__ElkSingleEdgeSection__Group_1_0_3__1__Impl"
    // InternalElkGraph.g:3082:1: rule__ElkSingleEdgeSection__Group_1_0_3__1__Impl : ( ':' ) ;
    public final void rule__ElkSingleEdgeSection__Group_1_0_3__1__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:3086:1: ( ( ':' ) )
            // InternalElkGraph.g:3087:1: ( ':' )
            {
            // InternalElkGraph.g:3087:1: ( ':' )
            // InternalElkGraph.g:3088:2: ':'
            {
             before(grammarAccess.getElkSingleEdgeSectionAccess().getColonKeyword_1_0_3_1()); 
            match(input,21,FOLLOW_2); 
             after(grammarAccess.getElkSingleEdgeSectionAccess().getColonKeyword_1_0_3_1()); 

            }


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ElkSingleEdgeSection__Group_1_0_3__1__Impl"


    // $ANTLR start "rule__ElkSingleEdgeSection__Group_1_0_3__2"
    // InternalElkGraph.g:3097:1: rule__ElkSingleEdgeSection__Group_1_0_3__2 : rule__ElkSingleEdgeSection__Group_1_0_3__2__Impl rule__ElkSingleEdgeSection__Group_1_0_3__3 ;
    public final void rule__ElkSingleEdgeSection__Group_1_0_3__2() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:3101:1: ( rule__ElkSingleEdgeSection__Group_1_0_3__2__Impl rule__ElkSingleEdgeSection__Group_1_0_3__3 )
            // InternalElkGraph.g:3102:2: rule__ElkSingleEdgeSection__Group_1_0_3__2__Impl rule__ElkSingleEdgeSection__Group_1_0_3__3
            {
            pushFollow(FOLLOW_18);
            rule__ElkSingleEdgeSection__Group_1_0_3__2__Impl();

            state._fsp--;

            pushFollow(FOLLOW_2);
            rule__ElkSingleEdgeSection__Group_1_0_3__3();

            state._fsp--;


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ElkSingleEdgeSection__Group_1_0_3__2"


    // $ANTLR start "rule__ElkSingleEdgeSection__Group_1_0_3__2__Impl"
    // InternalElkGraph.g:3109:1: rule__ElkSingleEdgeSection__Group_1_0_3__2__Impl : ( ( rule__ElkSingleEdgeSection__EndXAssignment_1_0_3_2 ) ) ;
    public final void rule__ElkSingleEdgeSection__Group_1_0_3__2__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:3113:1: ( ( ( rule__ElkSingleEdgeSection__EndXAssignment_1_0_3_2 ) ) )
            // InternalElkGraph.g:3114:1: ( ( rule__ElkSingleEdgeSection__EndXAssignment_1_0_3_2 ) )
            {
            // InternalElkGraph.g:3114:1: ( ( rule__ElkSingleEdgeSection__EndXAssignment_1_0_3_2 ) )
            // InternalElkGraph.g:3115:2: ( rule__ElkSingleEdgeSection__EndXAssignment_1_0_3_2 )
            {
             before(grammarAccess.getElkSingleEdgeSectionAccess().getEndXAssignment_1_0_3_2()); 
            // InternalElkGraph.g:3116:2: ( rule__ElkSingleEdgeSection__EndXAssignment_1_0_3_2 )
            // InternalElkGraph.g:3116:3: rule__ElkSingleEdgeSection__EndXAssignment_1_0_3_2
            {
            pushFollow(FOLLOW_2);
            rule__ElkSingleEdgeSection__EndXAssignment_1_0_3_2();

            state._fsp--;


            }

             after(grammarAccess.getElkSingleEdgeSectionAccess().getEndXAssignment_1_0_3_2()); 

            }


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ElkSingleEdgeSection__Group_1_0_3__2__Impl"


    // $ANTLR start "rule__ElkSingleEdgeSection__Group_1_0_3__3"
    // InternalElkGraph.g:3124:1: rule__ElkSingleEdgeSection__Group_1_0_3__3 : rule__ElkSingleEdgeSection__Group_1_0_3__3__Impl rule__ElkSingleEdgeSection__Group_1_0_3__4 ;
    public final void rule__ElkSingleEdgeSection__Group_1_0_3__3() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:3128:1: ( rule__ElkSingleEdgeSection__Group_1_0_3__3__Impl rule__ElkSingleEdgeSection__Group_1_0_3__4 )
            // InternalElkGraph.g:3129:2: rule__ElkSingleEdgeSection__Group_1_0_3__3__Impl rule__ElkSingleEdgeSection__Group_1_0_3__4
            {
            pushFollow(FOLLOW_17);
            rule__ElkSingleEdgeSection__Group_1_0_3__3__Impl();

            state._fsp--;

            pushFollow(FOLLOW_2);
            rule__ElkSingleEdgeSection__Group_1_0_3__4();

            state._fsp--;


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ElkSingleEdgeSection__Group_1_0_3__3"


    // $ANTLR start "rule__ElkSingleEdgeSection__Group_1_0_3__3__Impl"
    // InternalElkGraph.g:3136:1: rule__ElkSingleEdgeSection__Group_1_0_3__3__Impl : ( ',' ) ;
    public final void rule__ElkSingleEdgeSection__Group_1_0_3__3__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:3140:1: ( ( ',' ) )
            // InternalElkGraph.g:3141:1: ( ',' )
            {
            // InternalElkGraph.g:3141:1: ( ',' )
            // InternalElkGraph.g:3142:2: ','
            {
             before(grammarAccess.getElkSingleEdgeSectionAccess().getCommaKeyword_1_0_3_3()); 
            match(input,27,FOLLOW_2); 
             after(grammarAccess.getElkSingleEdgeSectionAccess().getCommaKeyword_1_0_3_3()); 

            }


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ElkSingleEdgeSection__Group_1_0_3__3__Impl"


    // $ANTLR start "rule__ElkSingleEdgeSection__Group_1_0_3__4"
    // InternalElkGraph.g:3151:1: rule__ElkSingleEdgeSection__Group_1_0_3__4 : rule__ElkSingleEdgeSection__Group_1_0_3__4__Impl ;
    public final void rule__ElkSingleEdgeSection__Group_1_0_3__4() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:3155:1: ( rule__ElkSingleEdgeSection__Group_1_0_3__4__Impl )
            // InternalElkGraph.g:3156:2: rule__ElkSingleEdgeSection__Group_1_0_3__4__Impl
            {
            pushFollow(FOLLOW_2);
            rule__ElkSingleEdgeSection__Group_1_0_3__4__Impl();

            state._fsp--;


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ElkSingleEdgeSection__Group_1_0_3__4"


    // $ANTLR start "rule__ElkSingleEdgeSection__Group_1_0_3__4__Impl"
    // InternalElkGraph.g:3162:1: rule__ElkSingleEdgeSection__Group_1_0_3__4__Impl : ( ( rule__ElkSingleEdgeSection__EndYAssignment_1_0_3_4 ) ) ;
    public final void rule__ElkSingleEdgeSection__Group_1_0_3__4__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:3166:1: ( ( ( rule__ElkSingleEdgeSection__EndYAssignment_1_0_3_4 ) ) )
            // InternalElkGraph.g:3167:1: ( ( rule__ElkSingleEdgeSection__EndYAssignment_1_0_3_4 ) )
            {
            // InternalElkGraph.g:3167:1: ( ( rule__ElkSingleEdgeSection__EndYAssignment_1_0_3_4 ) )
            // InternalElkGraph.g:3168:2: ( rule__ElkSingleEdgeSection__EndYAssignment_1_0_3_4 )
            {
             before(grammarAccess.getElkSingleEdgeSectionAccess().getEndYAssignment_1_0_3_4()); 
            // InternalElkGraph.g:3169:2: ( rule__ElkSingleEdgeSection__EndYAssignment_1_0_3_4 )
            // InternalElkGraph.g:3169:3: rule__ElkSingleEdgeSection__EndYAssignment_1_0_3_4
            {
            pushFollow(FOLLOW_2);
            rule__ElkSingleEdgeSection__EndYAssignment_1_0_3_4();

            state._fsp--;


            }

             after(grammarAccess.getElkSingleEdgeSectionAccess().getEndYAssignment_1_0_3_4()); 

            }


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ElkSingleEdgeSection__Group_1_0_3__4__Impl"


    // $ANTLR start "rule__ElkSingleEdgeSection__Group_1_1__0"
    // InternalElkGraph.g:3178:1: rule__ElkSingleEdgeSection__Group_1_1__0 : rule__ElkSingleEdgeSection__Group_1_1__0__Impl rule__ElkSingleEdgeSection__Group_1_1__1 ;
    public final void rule__ElkSingleEdgeSection__Group_1_1__0() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:3182:1: ( rule__ElkSingleEdgeSection__Group_1_1__0__Impl rule__ElkSingleEdgeSection__Group_1_1__1 )
            // InternalElkGraph.g:3183:2: rule__ElkSingleEdgeSection__Group_1_1__0__Impl rule__ElkSingleEdgeSection__Group_1_1__1
            {
            pushFollow(FOLLOW_11);
            rule__ElkSingleEdgeSection__Group_1_1__0__Impl();

            state._fsp--;

            pushFollow(FOLLOW_2);
            rule__ElkSingleEdgeSection__Group_1_1__1();

            state._fsp--;


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ElkSingleEdgeSection__Group_1_1__0"


    // $ANTLR start "rule__ElkSingleEdgeSection__Group_1_1__0__Impl"
    // InternalElkGraph.g:3190:1: rule__ElkSingleEdgeSection__Group_1_1__0__Impl : ( 'bends' ) ;
    public final void rule__ElkSingleEdgeSection__Group_1_1__0__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:3194:1: ( ( 'bends' ) )
            // InternalElkGraph.g:3195:1: ( 'bends' )
            {
            // InternalElkGraph.g:3195:1: ( 'bends' )
            // InternalElkGraph.g:3196:2: 'bends'
            {
             before(grammarAccess.getElkSingleEdgeSectionAccess().getBendsKeyword_1_1_0()); 
            match(input,35,FOLLOW_2); 
             after(grammarAccess.getElkSingleEdgeSectionAccess().getBendsKeyword_1_1_0()); 

            }


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ElkSingleEdgeSection__Group_1_1__0__Impl"


    // $ANTLR start "rule__ElkSingleEdgeSection__Group_1_1__1"
    // InternalElkGraph.g:3205:1: rule__ElkSingleEdgeSection__Group_1_1__1 : rule__ElkSingleEdgeSection__Group_1_1__1__Impl rule__ElkSingleEdgeSection__Group_1_1__2 ;
    public final void rule__ElkSingleEdgeSection__Group_1_1__1() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:3209:1: ( rule__ElkSingleEdgeSection__Group_1_1__1__Impl rule__ElkSingleEdgeSection__Group_1_1__2 )
            // InternalElkGraph.g:3210:2: rule__ElkSingleEdgeSection__Group_1_1__1__Impl rule__ElkSingleEdgeSection__Group_1_1__2
            {
            pushFollow(FOLLOW_17);
            rule__ElkSingleEdgeSection__Group_1_1__1__Impl();

            state._fsp--;

            pushFollow(FOLLOW_2);
            rule__ElkSingleEdgeSection__Group_1_1__2();

            state._fsp--;


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ElkSingleEdgeSection__Group_1_1__1"


    // $ANTLR start "rule__ElkSingleEdgeSection__Group_1_1__1__Impl"
    // InternalElkGraph.g:3217:1: rule__ElkSingleEdgeSection__Group_1_1__1__Impl : ( ':' ) ;
    public final void rule__ElkSingleEdgeSection__Group_1_1__1__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:3221:1: ( ( ':' ) )
            // InternalElkGraph.g:3222:1: ( ':' )
            {
            // InternalElkGraph.g:3222:1: ( ':' )
            // InternalElkGraph.g:3223:2: ':'
            {
             before(grammarAccess.getElkSingleEdgeSectionAccess().getColonKeyword_1_1_1()); 
            match(input,21,FOLLOW_2); 
             after(grammarAccess.getElkSingleEdgeSectionAccess().getColonKeyword_1_1_1()); 

            }


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ElkSingleEdgeSection__Group_1_1__1__Impl"


    // $ANTLR start "rule__ElkSingleEdgeSection__Group_1_1__2"
    // InternalElkGraph.g:3232:1: rule__ElkSingleEdgeSection__Group_1_1__2 : rule__ElkSingleEdgeSection__Group_1_1__2__Impl rule__ElkSingleEdgeSection__Group_1_1__3 ;
    public final void rule__ElkSingleEdgeSection__Group_1_1__2() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:3236:1: ( rule__ElkSingleEdgeSection__Group_1_1__2__Impl rule__ElkSingleEdgeSection__Group_1_1__3 )
            // InternalElkGraph.g:3237:2: rule__ElkSingleEdgeSection__Group_1_1__2__Impl rule__ElkSingleEdgeSection__Group_1_1__3
            {
            pushFollow(FOLLOW_25);
            rule__ElkSingleEdgeSection__Group_1_1__2__Impl();

            state._fsp--;

            pushFollow(FOLLOW_2);
            rule__ElkSingleEdgeSection__Group_1_1__3();

            state._fsp--;


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ElkSingleEdgeSection__Group_1_1__2"


    // $ANTLR start "rule__ElkSingleEdgeSection__Group_1_1__2__Impl"
    // InternalElkGraph.g:3244:1: rule__ElkSingleEdgeSection__Group_1_1__2__Impl : ( ( rule__ElkSingleEdgeSection__BendPointsAssignment_1_1_2 ) ) ;
    public final void rule__ElkSingleEdgeSection__Group_1_1__2__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:3248:1: ( ( ( rule__ElkSingleEdgeSection__BendPointsAssignment_1_1_2 ) ) )
            // InternalElkGraph.g:3249:1: ( ( rule__ElkSingleEdgeSection__BendPointsAssignment_1_1_2 ) )
            {
            // InternalElkGraph.g:3249:1: ( ( rule__ElkSingleEdgeSection__BendPointsAssignment_1_1_2 ) )
            // InternalElkGraph.g:3250:2: ( rule__ElkSingleEdgeSection__BendPointsAssignment_1_1_2 )
            {
             before(grammarAccess.getElkSingleEdgeSectionAccess().getBendPointsAssignment_1_1_2()); 
            // InternalElkGraph.g:3251:2: ( rule__ElkSingleEdgeSection__BendPointsAssignment_1_1_2 )
            // InternalElkGraph.g:3251:3: rule__ElkSingleEdgeSection__BendPointsAssignment_1_1_2
            {
            pushFollow(FOLLOW_2);
            rule__ElkSingleEdgeSection__BendPointsAssignment_1_1_2();

            state._fsp--;


            }

             after(grammarAccess.getElkSingleEdgeSectionAccess().getBendPointsAssignment_1_1_2()); 

            }


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ElkSingleEdgeSection__Group_1_1__2__Impl"


    // $ANTLR start "rule__ElkSingleEdgeSection__Group_1_1__3"
    // InternalElkGraph.g:3259:1: rule__ElkSingleEdgeSection__Group_1_1__3 : rule__ElkSingleEdgeSection__Group_1_1__3__Impl ;
    public final void rule__ElkSingleEdgeSection__Group_1_1__3() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:3263:1: ( rule__ElkSingleEdgeSection__Group_1_1__3__Impl )
            // InternalElkGraph.g:3264:2: rule__ElkSingleEdgeSection__Group_1_1__3__Impl
            {
            pushFollow(FOLLOW_2);
            rule__ElkSingleEdgeSection__Group_1_1__3__Impl();

            state._fsp--;


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ElkSingleEdgeSection__Group_1_1__3"


    // $ANTLR start "rule__ElkSingleEdgeSection__Group_1_1__3__Impl"
    // InternalElkGraph.g:3270:1: rule__ElkSingleEdgeSection__Group_1_1__3__Impl : ( ( rule__ElkSingleEdgeSection__Group_1_1_3__0 )* ) ;
    public final void rule__ElkSingleEdgeSection__Group_1_1__3__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:3274:1: ( ( ( rule__ElkSingleEdgeSection__Group_1_1_3__0 )* ) )
            // InternalElkGraph.g:3275:1: ( ( rule__ElkSingleEdgeSection__Group_1_1_3__0 )* )
            {
            // InternalElkGraph.g:3275:1: ( ( rule__ElkSingleEdgeSection__Group_1_1_3__0 )* )
            // InternalElkGraph.g:3276:2: ( rule__ElkSingleEdgeSection__Group_1_1_3__0 )*
            {
             before(grammarAccess.getElkSingleEdgeSectionAccess().getGroup_1_1_3()); 
            // InternalElkGraph.g:3277:2: ( rule__ElkSingleEdgeSection__Group_1_1_3__0 )*
            loop35:
            do {
                int alt35=2;
                int LA35_0 = input.LA(1);

                if ( (LA35_0==36) ) {
                    alt35=1;
                }


                switch (alt35) {
            	case 1 :
            	    // InternalElkGraph.g:3277:3: rule__ElkSingleEdgeSection__Group_1_1_3__0
            	    {
            	    pushFollow(FOLLOW_26);
            	    rule__ElkSingleEdgeSection__Group_1_1_3__0();

            	    state._fsp--;


            	    }
            	    break;

            	default :
            	    break loop35;
                }
            } while (true);

             after(grammarAccess.getElkSingleEdgeSectionAccess().getGroup_1_1_3()); 

            }


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ElkSingleEdgeSection__Group_1_1__3__Impl"


    // $ANTLR start "rule__ElkSingleEdgeSection__Group_1_1_3__0"
    // InternalElkGraph.g:3286:1: rule__ElkSingleEdgeSection__Group_1_1_3__0 : rule__ElkSingleEdgeSection__Group_1_1_3__0__Impl rule__ElkSingleEdgeSection__Group_1_1_3__1 ;
    public final void rule__ElkSingleEdgeSection__Group_1_1_3__0() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:3290:1: ( rule__ElkSingleEdgeSection__Group_1_1_3__0__Impl rule__ElkSingleEdgeSection__Group_1_1_3__1 )
            // InternalElkGraph.g:3291:2: rule__ElkSingleEdgeSection__Group_1_1_3__0__Impl rule__ElkSingleEdgeSection__Group_1_1_3__1
            {
            pushFollow(FOLLOW_17);
            rule__ElkSingleEdgeSection__Group_1_1_3__0__Impl();

            state._fsp--;

            pushFollow(FOLLOW_2);
            rule__ElkSingleEdgeSection__Group_1_1_3__1();

            state._fsp--;


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ElkSingleEdgeSection__Group_1_1_3__0"


    // $ANTLR start "rule__ElkSingleEdgeSection__Group_1_1_3__0__Impl"
    // InternalElkGraph.g:3298:1: rule__ElkSingleEdgeSection__Group_1_1_3__0__Impl : ( '|' ) ;
    public final void rule__ElkSingleEdgeSection__Group_1_1_3__0__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:3302:1: ( ( '|' ) )
            // InternalElkGraph.g:3303:1: ( '|' )
            {
            // InternalElkGraph.g:3303:1: ( '|' )
            // InternalElkGraph.g:3304:2: '|'
            {
             before(grammarAccess.getElkSingleEdgeSectionAccess().getVerticalLineKeyword_1_1_3_0()); 
            match(input,36,FOLLOW_2); 
             after(grammarAccess.getElkSingleEdgeSectionAccess().getVerticalLineKeyword_1_1_3_0()); 

            }


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ElkSingleEdgeSection__Group_1_1_3__0__Impl"


    // $ANTLR start "rule__ElkSingleEdgeSection__Group_1_1_3__1"
    // InternalElkGraph.g:3313:1: rule__ElkSingleEdgeSection__Group_1_1_3__1 : rule__ElkSingleEdgeSection__Group_1_1_3__1__Impl ;
    public final void rule__ElkSingleEdgeSection__Group_1_1_3__1() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:3317:1: ( rule__ElkSingleEdgeSection__Group_1_1_3__1__Impl )
            // InternalElkGraph.g:3318:2: rule__ElkSingleEdgeSection__Group_1_1_3__1__Impl
            {
            pushFollow(FOLLOW_2);
            rule__ElkSingleEdgeSection__Group_1_1_3__1__Impl();

            state._fsp--;


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ElkSingleEdgeSection__Group_1_1_3__1"


    // $ANTLR start "rule__ElkSingleEdgeSection__Group_1_1_3__1__Impl"
    // InternalElkGraph.g:3324:1: rule__ElkSingleEdgeSection__Group_1_1_3__1__Impl : ( ( rule__ElkSingleEdgeSection__BendPointsAssignment_1_1_3_1 ) ) ;
    public final void rule__ElkSingleEdgeSection__Group_1_1_3__1__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:3328:1: ( ( ( rule__ElkSingleEdgeSection__BendPointsAssignment_1_1_3_1 ) ) )
            // InternalElkGraph.g:3329:1: ( ( rule__ElkSingleEdgeSection__BendPointsAssignment_1_1_3_1 ) )
            {
            // InternalElkGraph.g:3329:1: ( ( rule__ElkSingleEdgeSection__BendPointsAssignment_1_1_3_1 ) )
            // InternalElkGraph.g:3330:2: ( rule__ElkSingleEdgeSection__BendPointsAssignment_1_1_3_1 )
            {
             before(grammarAccess.getElkSingleEdgeSectionAccess().getBendPointsAssignment_1_1_3_1()); 
            // InternalElkGraph.g:3331:2: ( rule__ElkSingleEdgeSection__BendPointsAssignment_1_1_3_1 )
            // InternalElkGraph.g:3331:3: rule__ElkSingleEdgeSection__BendPointsAssignment_1_1_3_1
            {
            pushFollow(FOLLOW_2);
            rule__ElkSingleEdgeSection__BendPointsAssignment_1_1_3_1();

            state._fsp--;


            }

             after(grammarAccess.getElkSingleEdgeSectionAccess().getBendPointsAssignment_1_1_3_1()); 

            }


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ElkSingleEdgeSection__Group_1_1_3__1__Impl"


    // $ANTLR start "rule__ElkEdgeSection__Group__0"
    // InternalElkGraph.g:3340:1: rule__ElkEdgeSection__Group__0 : rule__ElkEdgeSection__Group__0__Impl rule__ElkEdgeSection__Group__1 ;
    public final void rule__ElkEdgeSection__Group__0() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:3344:1: ( rule__ElkEdgeSection__Group__0__Impl rule__ElkEdgeSection__Group__1 )
            // InternalElkGraph.g:3345:2: rule__ElkEdgeSection__Group__0__Impl rule__ElkEdgeSection__Group__1
            {
            pushFollow(FOLLOW_7);
            rule__ElkEdgeSection__Group__0__Impl();

            state._fsp--;

            pushFollow(FOLLOW_2);
            rule__ElkEdgeSection__Group__1();

            state._fsp--;


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ElkEdgeSection__Group__0"


    // $ANTLR start "rule__ElkEdgeSection__Group__0__Impl"
    // InternalElkGraph.g:3352:1: rule__ElkEdgeSection__Group__0__Impl : ( 'section' ) ;
    public final void rule__ElkEdgeSection__Group__0__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:3356:1: ( ( 'section' ) )
            // InternalElkGraph.g:3357:1: ( 'section' )
            {
            // InternalElkGraph.g:3357:1: ( 'section' )
            // InternalElkGraph.g:3358:2: 'section'
            {
             before(grammarAccess.getElkEdgeSectionAccess().getSectionKeyword_0()); 
            match(input,37,FOLLOW_2); 
             after(grammarAccess.getElkEdgeSectionAccess().getSectionKeyword_0()); 

            }


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ElkEdgeSection__Group__0__Impl"


    // $ANTLR start "rule__ElkEdgeSection__Group__1"
    // InternalElkGraph.g:3367:1: rule__ElkEdgeSection__Group__1 : rule__ElkEdgeSection__Group__1__Impl rule__ElkEdgeSection__Group__2 ;
    public final void rule__ElkEdgeSection__Group__1() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:3371:1: ( rule__ElkEdgeSection__Group__1__Impl rule__ElkEdgeSection__Group__2 )
            // InternalElkGraph.g:3372:2: rule__ElkEdgeSection__Group__1__Impl rule__ElkEdgeSection__Group__2
            {
            pushFollow(FOLLOW_27);
            rule__ElkEdgeSection__Group__1__Impl();

            state._fsp--;

            pushFollow(FOLLOW_2);
            rule__ElkEdgeSection__Group__2();

            state._fsp--;


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ElkEdgeSection__Group__1"


    // $ANTLR start "rule__ElkEdgeSection__Group__1__Impl"
    // InternalElkGraph.g:3379:1: rule__ElkEdgeSection__Group__1__Impl : ( ( rule__ElkEdgeSection__IdentifierAssignment_1 ) ) ;
    public final void rule__ElkEdgeSection__Group__1__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:3383:1: ( ( ( rule__ElkEdgeSection__IdentifierAssignment_1 ) ) )
            // InternalElkGraph.g:3384:1: ( ( rule__ElkEdgeSection__IdentifierAssignment_1 ) )
            {
            // InternalElkGraph.g:3384:1: ( ( rule__ElkEdgeSection__IdentifierAssignment_1 ) )
            // InternalElkGraph.g:3385:2: ( rule__ElkEdgeSection__IdentifierAssignment_1 )
            {
             before(grammarAccess.getElkEdgeSectionAccess().getIdentifierAssignment_1()); 
            // InternalElkGraph.g:3386:2: ( rule__ElkEdgeSection__IdentifierAssignment_1 )
            // InternalElkGraph.g:3386:3: rule__ElkEdgeSection__IdentifierAssignment_1
            {
            pushFollow(FOLLOW_2);
            rule__ElkEdgeSection__IdentifierAssignment_1();

            state._fsp--;


            }

             after(grammarAccess.getElkEdgeSectionAccess().getIdentifierAssignment_1()); 

            }


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ElkEdgeSection__Group__1__Impl"


    // $ANTLR start "rule__ElkEdgeSection__Group__2"
    // InternalElkGraph.g:3394:1: rule__ElkEdgeSection__Group__2 : rule__ElkEdgeSection__Group__2__Impl rule__ElkEdgeSection__Group__3 ;
    public final void rule__ElkEdgeSection__Group__2() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:3398:1: ( rule__ElkEdgeSection__Group__2__Impl rule__ElkEdgeSection__Group__3 )
            // InternalElkGraph.g:3399:2: rule__ElkEdgeSection__Group__2__Impl rule__ElkEdgeSection__Group__3
            {
            pushFollow(FOLLOW_27);
            rule__ElkEdgeSection__Group__2__Impl();

            state._fsp--;

            pushFollow(FOLLOW_2);
            rule__ElkEdgeSection__Group__3();

            state._fsp--;


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ElkEdgeSection__Group__2"


    // $ANTLR start "rule__ElkEdgeSection__Group__2__Impl"
    // InternalElkGraph.g:3406:1: rule__ElkEdgeSection__Group__2__Impl : ( ( rule__ElkEdgeSection__Group_2__0 )? ) ;
    public final void rule__ElkEdgeSection__Group__2__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:3410:1: ( ( ( rule__ElkEdgeSection__Group_2__0 )? ) )
            // InternalElkGraph.g:3411:1: ( ( rule__ElkEdgeSection__Group_2__0 )? )
            {
            // InternalElkGraph.g:3411:1: ( ( rule__ElkEdgeSection__Group_2__0 )? )
            // InternalElkGraph.g:3412:2: ( rule__ElkEdgeSection__Group_2__0 )?
            {
             before(grammarAccess.getElkEdgeSectionAccess().getGroup_2()); 
            // InternalElkGraph.g:3413:2: ( rule__ElkEdgeSection__Group_2__0 )?
            int alt36=2;
            int LA36_0 = input.LA(1);

            if ( (LA36_0==30) ) {
                alt36=1;
            }
            switch (alt36) {
                case 1 :
                    // InternalElkGraph.g:3413:3: rule__ElkEdgeSection__Group_2__0
                    {
                    pushFollow(FOLLOW_2);
                    rule__ElkEdgeSection__Group_2__0();

                    state._fsp--;


                    }
                    break;

            }

             after(grammarAccess.getElkEdgeSectionAccess().getGroup_2()); 

            }


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ElkEdgeSection__Group__2__Impl"


    // $ANTLR start "rule__ElkEdgeSection__Group__3"
    // InternalElkGraph.g:3421:1: rule__ElkEdgeSection__Group__3 : rule__ElkEdgeSection__Group__3__Impl rule__ElkEdgeSection__Group__4 ;
    public final void rule__ElkEdgeSection__Group__3() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:3425:1: ( rule__ElkEdgeSection__Group__3__Impl rule__ElkEdgeSection__Group__4 )
            // InternalElkGraph.g:3426:2: rule__ElkEdgeSection__Group__3__Impl rule__ElkEdgeSection__Group__4
            {
            pushFollow(FOLLOW_23);
            rule__ElkEdgeSection__Group__3__Impl();

            state._fsp--;

            pushFollow(FOLLOW_2);
            rule__ElkEdgeSection__Group__4();

            state._fsp--;


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ElkEdgeSection__Group__3"


    // $ANTLR start "rule__ElkEdgeSection__Group__3__Impl"
    // InternalElkGraph.g:3433:1: rule__ElkEdgeSection__Group__3__Impl : ( '[' ) ;
    public final void rule__ElkEdgeSection__Group__3__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:3437:1: ( ( '[' ) )
            // InternalElkGraph.g:3438:1: ( '[' )
            {
            // InternalElkGraph.g:3438:1: ( '[' )
            // InternalElkGraph.g:3439:2: '['
            {
             before(grammarAccess.getElkEdgeSectionAccess().getLeftSquareBracketKeyword_3()); 
            match(input,24,FOLLOW_2); 
             after(grammarAccess.getElkEdgeSectionAccess().getLeftSquareBracketKeyword_3()); 

            }


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ElkEdgeSection__Group__3__Impl"


    // $ANTLR start "rule__ElkEdgeSection__Group__4"
    // InternalElkGraph.g:3448:1: rule__ElkEdgeSection__Group__4 : rule__ElkEdgeSection__Group__4__Impl rule__ElkEdgeSection__Group__5 ;
    public final void rule__ElkEdgeSection__Group__4() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:3452:1: ( rule__ElkEdgeSection__Group__4__Impl rule__ElkEdgeSection__Group__5 )
            // InternalElkGraph.g:3453:2: rule__ElkEdgeSection__Group__4__Impl rule__ElkEdgeSection__Group__5
            {
            pushFollow(FOLLOW_16);
            rule__ElkEdgeSection__Group__4__Impl();

            state._fsp--;

            pushFollow(FOLLOW_2);
            rule__ElkEdgeSection__Group__5();

            state._fsp--;


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ElkEdgeSection__Group__4"


    // $ANTLR start "rule__ElkEdgeSection__Group__4__Impl"
    // InternalElkGraph.g:3460:1: rule__ElkEdgeSection__Group__4__Impl : ( ( rule__ElkEdgeSection__Group_4__0 ) ) ;
    public final void rule__ElkEdgeSection__Group__4__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:3464:1: ( ( ( rule__ElkEdgeSection__Group_4__0 ) ) )
            // InternalElkGraph.g:3465:1: ( ( rule__ElkEdgeSection__Group_4__0 ) )
            {
            // InternalElkGraph.g:3465:1: ( ( rule__ElkEdgeSection__Group_4__0 ) )
            // InternalElkGraph.g:3466:2: ( rule__ElkEdgeSection__Group_4__0 )
            {
             before(grammarAccess.getElkEdgeSectionAccess().getGroup_4()); 
            // InternalElkGraph.g:3467:2: ( rule__ElkEdgeSection__Group_4__0 )
            // InternalElkGraph.g:3467:3: rule__ElkEdgeSection__Group_4__0
            {
            pushFollow(FOLLOW_2);
            rule__ElkEdgeSection__Group_4__0();

            state._fsp--;


            }

             after(grammarAccess.getElkEdgeSectionAccess().getGroup_4()); 

            }


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ElkEdgeSection__Group__4__Impl"


    // $ANTLR start "rule__ElkEdgeSection__Group__5"
    // InternalElkGraph.g:3475:1: rule__ElkEdgeSection__Group__5 : rule__ElkEdgeSection__Group__5__Impl ;
    public final void rule__ElkEdgeSection__Group__5() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:3479:1: ( rule__ElkEdgeSection__Group__5__Impl )
            // InternalElkGraph.g:3480:2: rule__ElkEdgeSection__Group__5__Impl
            {
            pushFollow(FOLLOW_2);
            rule__ElkEdgeSection__Group__5__Impl();

            state._fsp--;


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ElkEdgeSection__Group__5"


    // $ANTLR start "rule__ElkEdgeSection__Group__5__Impl"
    // InternalElkGraph.g:3486:1: rule__ElkEdgeSection__Group__5__Impl : ( ']' ) ;
    public final void rule__ElkEdgeSection__Group__5__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:3490:1: ( ( ']' ) )
            // InternalElkGraph.g:3491:1: ( ']' )
            {
            // InternalElkGraph.g:3491:1: ( ']' )
            // InternalElkGraph.g:3492:2: ']'
            {
             before(grammarAccess.getElkEdgeSectionAccess().getRightSquareBracketKeyword_5()); 
            match(input,25,FOLLOW_2); 
             after(grammarAccess.getElkEdgeSectionAccess().getRightSquareBracketKeyword_5()); 

            }


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ElkEdgeSection__Group__5__Impl"


    // $ANTLR start "rule__ElkEdgeSection__Group_2__0"
    // InternalElkGraph.g:3502:1: rule__ElkEdgeSection__Group_2__0 : rule__ElkEdgeSection__Group_2__0__Impl rule__ElkEdgeSection__Group_2__1 ;
    public final void rule__ElkEdgeSection__Group_2__0() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:3506:1: ( rule__ElkEdgeSection__Group_2__0__Impl rule__ElkEdgeSection__Group_2__1 )
            // InternalElkGraph.g:3507:2: rule__ElkEdgeSection__Group_2__0__Impl rule__ElkEdgeSection__Group_2__1
            {
            pushFollow(FOLLOW_7);
            rule__ElkEdgeSection__Group_2__0__Impl();

            state._fsp--;

            pushFollow(FOLLOW_2);
            rule__ElkEdgeSection__Group_2__1();

            state._fsp--;


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ElkEdgeSection__Group_2__0"


    // $ANTLR start "rule__ElkEdgeSection__Group_2__0__Impl"
    // InternalElkGraph.g:3514:1: rule__ElkEdgeSection__Group_2__0__Impl : ( '->' ) ;
    public final void rule__ElkEdgeSection__Group_2__0__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:3518:1: ( ( '->' ) )
            // InternalElkGraph.g:3519:1: ( '->' )
            {
            // InternalElkGraph.g:3519:1: ( '->' )
            // InternalElkGraph.g:3520:2: '->'
            {
             before(grammarAccess.getElkEdgeSectionAccess().getHyphenMinusGreaterThanSignKeyword_2_0()); 
            match(input,30,FOLLOW_2); 
             after(grammarAccess.getElkEdgeSectionAccess().getHyphenMinusGreaterThanSignKeyword_2_0()); 

            }


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ElkEdgeSection__Group_2__0__Impl"


    // $ANTLR start "rule__ElkEdgeSection__Group_2__1"
    // InternalElkGraph.g:3529:1: rule__ElkEdgeSection__Group_2__1 : rule__ElkEdgeSection__Group_2__1__Impl rule__ElkEdgeSection__Group_2__2 ;
    public final void rule__ElkEdgeSection__Group_2__1() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:3533:1: ( rule__ElkEdgeSection__Group_2__1__Impl rule__ElkEdgeSection__Group_2__2 )
            // InternalElkGraph.g:3534:2: rule__ElkEdgeSection__Group_2__1__Impl rule__ElkEdgeSection__Group_2__2
            {
            pushFollow(FOLLOW_18);
            rule__ElkEdgeSection__Group_2__1__Impl();

            state._fsp--;

            pushFollow(FOLLOW_2);
            rule__ElkEdgeSection__Group_2__2();

            state._fsp--;


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ElkEdgeSection__Group_2__1"


    // $ANTLR start "rule__ElkEdgeSection__Group_2__1__Impl"
    // InternalElkGraph.g:3541:1: rule__ElkEdgeSection__Group_2__1__Impl : ( ( rule__ElkEdgeSection__OutgoingSectionsAssignment_2_1 ) ) ;
    public final void rule__ElkEdgeSection__Group_2__1__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:3545:1: ( ( ( rule__ElkEdgeSection__OutgoingSectionsAssignment_2_1 ) ) )
            // InternalElkGraph.g:3546:1: ( ( rule__ElkEdgeSection__OutgoingSectionsAssignment_2_1 ) )
            {
            // InternalElkGraph.g:3546:1: ( ( rule__ElkEdgeSection__OutgoingSectionsAssignment_2_1 ) )
            // InternalElkGraph.g:3547:2: ( rule__ElkEdgeSection__OutgoingSectionsAssignment_2_1 )
            {
             before(grammarAccess.getElkEdgeSectionAccess().getOutgoingSectionsAssignment_2_1()); 
            // InternalElkGraph.g:3548:2: ( rule__ElkEdgeSection__OutgoingSectionsAssignment_2_1 )
            // InternalElkGraph.g:3548:3: rule__ElkEdgeSection__OutgoingSectionsAssignment_2_1
            {
            pushFollow(FOLLOW_2);
            rule__ElkEdgeSection__OutgoingSectionsAssignment_2_1();

            state._fsp--;


            }

             after(grammarAccess.getElkEdgeSectionAccess().getOutgoingSectionsAssignment_2_1()); 

            }


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ElkEdgeSection__Group_2__1__Impl"


    // $ANTLR start "rule__ElkEdgeSection__Group_2__2"
    // InternalElkGraph.g:3556:1: rule__ElkEdgeSection__Group_2__2 : rule__ElkEdgeSection__Group_2__2__Impl ;
    public final void rule__ElkEdgeSection__Group_2__2() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:3560:1: ( rule__ElkEdgeSection__Group_2__2__Impl )
            // InternalElkGraph.g:3561:2: rule__ElkEdgeSection__Group_2__2__Impl
            {
            pushFollow(FOLLOW_2);
            rule__ElkEdgeSection__Group_2__2__Impl();

            state._fsp--;


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ElkEdgeSection__Group_2__2"


    // $ANTLR start "rule__ElkEdgeSection__Group_2__2__Impl"
    // InternalElkGraph.g:3567:1: rule__ElkEdgeSection__Group_2__2__Impl : ( ( rule__ElkEdgeSection__Group_2_2__0 )* ) ;
    public final void rule__ElkEdgeSection__Group_2__2__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:3571:1: ( ( ( rule__ElkEdgeSection__Group_2_2__0 )* ) )
            // InternalElkGraph.g:3572:1: ( ( rule__ElkEdgeSection__Group_2_2__0 )* )
            {
            // InternalElkGraph.g:3572:1: ( ( rule__ElkEdgeSection__Group_2_2__0 )* )
            // InternalElkGraph.g:3573:2: ( rule__ElkEdgeSection__Group_2_2__0 )*
            {
             before(grammarAccess.getElkEdgeSectionAccess().getGroup_2_2()); 
            // InternalElkGraph.g:3574:2: ( rule__ElkEdgeSection__Group_2_2__0 )*
            loop37:
            do {
                int alt37=2;
                int LA37_0 = input.LA(1);

                if ( (LA37_0==27) ) {
                    alt37=1;
                }


                switch (alt37) {
            	case 1 :
            	    // InternalElkGraph.g:3574:3: rule__ElkEdgeSection__Group_2_2__0
            	    {
            	    pushFollow(FOLLOW_20);
            	    rule__ElkEdgeSection__Group_2_2__0();

            	    state._fsp--;


            	    }
            	    break;

            	default :
            	    break loop37;
                }
            } while (true);

             after(grammarAccess.getElkEdgeSectionAccess().getGroup_2_2()); 

            }


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ElkEdgeSection__Group_2__2__Impl"


    // $ANTLR start "rule__ElkEdgeSection__Group_2_2__0"
    // InternalElkGraph.g:3583:1: rule__ElkEdgeSection__Group_2_2__0 : rule__ElkEdgeSection__Group_2_2__0__Impl rule__ElkEdgeSection__Group_2_2__1 ;
    public final void rule__ElkEdgeSection__Group_2_2__0() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:3587:1: ( rule__ElkEdgeSection__Group_2_2__0__Impl rule__ElkEdgeSection__Group_2_2__1 )
            // InternalElkGraph.g:3588:2: rule__ElkEdgeSection__Group_2_2__0__Impl rule__ElkEdgeSection__Group_2_2__1
            {
            pushFollow(FOLLOW_7);
            rule__ElkEdgeSection__Group_2_2__0__Impl();

            state._fsp--;

            pushFollow(FOLLOW_2);
            rule__ElkEdgeSection__Group_2_2__1();

            state._fsp--;


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ElkEdgeSection__Group_2_2__0"


    // $ANTLR start "rule__ElkEdgeSection__Group_2_2__0__Impl"
    // InternalElkGraph.g:3595:1: rule__ElkEdgeSection__Group_2_2__0__Impl : ( ',' ) ;
    public final void rule__ElkEdgeSection__Group_2_2__0__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:3599:1: ( ( ',' ) )
            // InternalElkGraph.g:3600:1: ( ',' )
            {
            // InternalElkGraph.g:3600:1: ( ',' )
            // InternalElkGraph.g:3601:2: ','
            {
             before(grammarAccess.getElkEdgeSectionAccess().getCommaKeyword_2_2_0()); 
            match(input,27,FOLLOW_2); 
             after(grammarAccess.getElkEdgeSectionAccess().getCommaKeyword_2_2_0()); 

            }


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ElkEdgeSection__Group_2_2__0__Impl"


    // $ANTLR start "rule__ElkEdgeSection__Group_2_2__1"
    // InternalElkGraph.g:3610:1: rule__ElkEdgeSection__Group_2_2__1 : rule__ElkEdgeSection__Group_2_2__1__Impl ;
    public final void rule__ElkEdgeSection__Group_2_2__1() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:3614:1: ( rule__ElkEdgeSection__Group_2_2__1__Impl )
            // InternalElkGraph.g:3615:2: rule__ElkEdgeSection__Group_2_2__1__Impl
            {
            pushFollow(FOLLOW_2);
            rule__ElkEdgeSection__Group_2_2__1__Impl();

            state._fsp--;


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ElkEdgeSection__Group_2_2__1"


    // $ANTLR start "rule__ElkEdgeSection__Group_2_2__1__Impl"
    // InternalElkGraph.g:3621:1: rule__ElkEdgeSection__Group_2_2__1__Impl : ( ( rule__ElkEdgeSection__OutgoingSectionsAssignment_2_2_1 ) ) ;
    public final void rule__ElkEdgeSection__Group_2_2__1__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:3625:1: ( ( ( rule__ElkEdgeSection__OutgoingSectionsAssignment_2_2_1 ) ) )
            // InternalElkGraph.g:3626:1: ( ( rule__ElkEdgeSection__OutgoingSectionsAssignment_2_2_1 ) )
            {
            // InternalElkGraph.g:3626:1: ( ( rule__ElkEdgeSection__OutgoingSectionsAssignment_2_2_1 ) )
            // InternalElkGraph.g:3627:2: ( rule__ElkEdgeSection__OutgoingSectionsAssignment_2_2_1 )
            {
             before(grammarAccess.getElkEdgeSectionAccess().getOutgoingSectionsAssignment_2_2_1()); 
            // InternalElkGraph.g:3628:2: ( rule__ElkEdgeSection__OutgoingSectionsAssignment_2_2_1 )
            // InternalElkGraph.g:3628:3: rule__ElkEdgeSection__OutgoingSectionsAssignment_2_2_1
            {
            pushFollow(FOLLOW_2);
            rule__ElkEdgeSection__OutgoingSectionsAssignment_2_2_1();

            state._fsp--;


            }

             after(grammarAccess.getElkEdgeSectionAccess().getOutgoingSectionsAssignment_2_2_1()); 

            }


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ElkEdgeSection__Group_2_2__1__Impl"


    // $ANTLR start "rule__ElkEdgeSection__Group_4__0"
    // InternalElkGraph.g:3637:1: rule__ElkEdgeSection__Group_4__0 : rule__ElkEdgeSection__Group_4__0__Impl rule__ElkEdgeSection__Group_4__1 ;
    public final void rule__ElkEdgeSection__Group_4__0() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:3641:1: ( rule__ElkEdgeSection__Group_4__0__Impl rule__ElkEdgeSection__Group_4__1 )
            // InternalElkGraph.g:3642:2: rule__ElkEdgeSection__Group_4__0__Impl rule__ElkEdgeSection__Group_4__1
            {
            pushFollow(FOLLOW_24);
            rule__ElkEdgeSection__Group_4__0__Impl();

            state._fsp--;

            pushFollow(FOLLOW_2);
            rule__ElkEdgeSection__Group_4__1();

            state._fsp--;


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ElkEdgeSection__Group_4__0"


    // $ANTLR start "rule__ElkEdgeSection__Group_4__0__Impl"
    // InternalElkGraph.g:3649:1: rule__ElkEdgeSection__Group_4__0__Impl : ( ( rule__ElkEdgeSection__UnorderedGroup_4_0 ) ) ;
    public final void rule__ElkEdgeSection__Group_4__0__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:3653:1: ( ( ( rule__ElkEdgeSection__UnorderedGroup_4_0 ) ) )
            // InternalElkGraph.g:3654:1: ( ( rule__ElkEdgeSection__UnorderedGroup_4_0 ) )
            {
            // InternalElkGraph.g:3654:1: ( ( rule__ElkEdgeSection__UnorderedGroup_4_0 ) )
            // InternalElkGraph.g:3655:2: ( rule__ElkEdgeSection__UnorderedGroup_4_0 )
            {
             before(grammarAccess.getElkEdgeSectionAccess().getUnorderedGroup_4_0()); 
            // InternalElkGraph.g:3656:2: ( rule__ElkEdgeSection__UnorderedGroup_4_0 )
            // InternalElkGraph.g:3656:3: rule__ElkEdgeSection__UnorderedGroup_4_0
            {
            pushFollow(FOLLOW_2);
            rule__ElkEdgeSection__UnorderedGroup_4_0();

            state._fsp--;


            }

             after(grammarAccess.getElkEdgeSectionAccess().getUnorderedGroup_4_0()); 

            }


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ElkEdgeSection__Group_4__0__Impl"


    // $ANTLR start "rule__ElkEdgeSection__Group_4__1"
    // InternalElkGraph.g:3664:1: rule__ElkEdgeSection__Group_4__1 : rule__ElkEdgeSection__Group_4__1__Impl rule__ElkEdgeSection__Group_4__2 ;
    public final void rule__ElkEdgeSection__Group_4__1() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:3668:1: ( rule__ElkEdgeSection__Group_4__1__Impl rule__ElkEdgeSection__Group_4__2 )
            // InternalElkGraph.g:3669:2: rule__ElkEdgeSection__Group_4__1__Impl rule__ElkEdgeSection__Group_4__2
            {
            pushFollow(FOLLOW_24);
            rule__ElkEdgeSection__Group_4__1__Impl();

            state._fsp--;

            pushFollow(FOLLOW_2);
            rule__ElkEdgeSection__Group_4__2();

            state._fsp--;


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ElkEdgeSection__Group_4__1"


    // $ANTLR start "rule__ElkEdgeSection__Group_4__1__Impl"
    // InternalElkGraph.g:3676:1: rule__ElkEdgeSection__Group_4__1__Impl : ( ( rule__ElkEdgeSection__Group_4_1__0 )? ) ;
    public final void rule__ElkEdgeSection__Group_4__1__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:3680:1: ( ( ( rule__ElkEdgeSection__Group_4_1__0 )? ) )
            // InternalElkGraph.g:3681:1: ( ( rule__ElkEdgeSection__Group_4_1__0 )? )
            {
            // InternalElkGraph.g:3681:1: ( ( rule__ElkEdgeSection__Group_4_1__0 )? )
            // InternalElkGraph.g:3682:2: ( rule__ElkEdgeSection__Group_4_1__0 )?
            {
             before(grammarAccess.getElkEdgeSectionAccess().getGroup_4_1()); 
            // InternalElkGraph.g:3683:2: ( rule__ElkEdgeSection__Group_4_1__0 )?
            int alt38=2;
            int LA38_0 = input.LA(1);

            if ( (LA38_0==35) ) {
                alt38=1;
            }
            switch (alt38) {
                case 1 :
                    // InternalElkGraph.g:3683:3: rule__ElkEdgeSection__Group_4_1__0
                    {
                    pushFollow(FOLLOW_2);
                    rule__ElkEdgeSection__Group_4_1__0();

                    state._fsp--;


                    }
                    break;

            }

             after(grammarAccess.getElkEdgeSectionAccess().getGroup_4_1()); 

            }


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ElkEdgeSection__Group_4__1__Impl"


    // $ANTLR start "rule__ElkEdgeSection__Group_4__2"
    // InternalElkGraph.g:3691:1: rule__ElkEdgeSection__Group_4__2 : rule__ElkEdgeSection__Group_4__2__Impl ;
    public final void rule__ElkEdgeSection__Group_4__2() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:3695:1: ( rule__ElkEdgeSection__Group_4__2__Impl )
            // InternalElkGraph.g:3696:2: rule__ElkEdgeSection__Group_4__2__Impl
            {
            pushFollow(FOLLOW_2);
            rule__ElkEdgeSection__Group_4__2__Impl();

            state._fsp--;


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ElkEdgeSection__Group_4__2"


    // $ANTLR start "rule__ElkEdgeSection__Group_4__2__Impl"
    // InternalElkGraph.g:3702:1: rule__ElkEdgeSection__Group_4__2__Impl : ( ( rule__ElkEdgeSection__PropertiesAssignment_4_2 )* ) ;
    public final void rule__ElkEdgeSection__Group_4__2__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:3706:1: ( ( ( rule__ElkEdgeSection__PropertiesAssignment_4_2 )* ) )
            // InternalElkGraph.g:3707:1: ( ( rule__ElkEdgeSection__PropertiesAssignment_4_2 )* )
            {
            // InternalElkGraph.g:3707:1: ( ( rule__ElkEdgeSection__PropertiesAssignment_4_2 )* )
            // InternalElkGraph.g:3708:2: ( rule__ElkEdgeSection__PropertiesAssignment_4_2 )*
            {
             before(grammarAccess.getElkEdgeSectionAccess().getPropertiesAssignment_4_2()); 
            // InternalElkGraph.g:3709:2: ( rule__ElkEdgeSection__PropertiesAssignment_4_2 )*
            loop39:
            do {
                int alt39=2;
                int LA39_0 = input.LA(1);

                if ( (LA39_0==RULE_ID) ) {
                    alt39=1;
                }


                switch (alt39) {
            	case 1 :
            	    // InternalElkGraph.g:3709:3: rule__ElkEdgeSection__PropertiesAssignment_4_2
            	    {
            	    pushFollow(FOLLOW_5);
            	    rule__ElkEdgeSection__PropertiesAssignment_4_2();

            	    state._fsp--;


            	    }
            	    break;

            	default :
            	    break loop39;
                }
            } while (true);

             after(grammarAccess.getElkEdgeSectionAccess().getPropertiesAssignment_4_2()); 

            }


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ElkEdgeSection__Group_4__2__Impl"


    // $ANTLR start "rule__ElkEdgeSection__Group_4_0_0__0"
    // InternalElkGraph.g:3718:1: rule__ElkEdgeSection__Group_4_0_0__0 : rule__ElkEdgeSection__Group_4_0_0__0__Impl rule__ElkEdgeSection__Group_4_0_0__1 ;
    public final void rule__ElkEdgeSection__Group_4_0_0__0() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:3722:1: ( rule__ElkEdgeSection__Group_4_0_0__0__Impl rule__ElkEdgeSection__Group_4_0_0__1 )
            // InternalElkGraph.g:3723:2: rule__ElkEdgeSection__Group_4_0_0__0__Impl rule__ElkEdgeSection__Group_4_0_0__1
            {
            pushFollow(FOLLOW_11);
            rule__ElkEdgeSection__Group_4_0_0__0__Impl();

            state._fsp--;

            pushFollow(FOLLOW_2);
            rule__ElkEdgeSection__Group_4_0_0__1();

            state._fsp--;


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ElkEdgeSection__Group_4_0_0__0"


    // $ANTLR start "rule__ElkEdgeSection__Group_4_0_0__0__Impl"
    // InternalElkGraph.g:3730:1: rule__ElkEdgeSection__Group_4_0_0__0__Impl : ( 'incoming' ) ;
    public final void rule__ElkEdgeSection__Group_4_0_0__0__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:3734:1: ( ( 'incoming' ) )
            // InternalElkGraph.g:3735:1: ( 'incoming' )
            {
            // InternalElkGraph.g:3735:1: ( 'incoming' )
            // InternalElkGraph.g:3736:2: 'incoming'
            {
             before(grammarAccess.getElkEdgeSectionAccess().getIncomingKeyword_4_0_0_0()); 
            match(input,31,FOLLOW_2); 
             after(grammarAccess.getElkEdgeSectionAccess().getIncomingKeyword_4_0_0_0()); 

            }


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ElkEdgeSection__Group_4_0_0__0__Impl"


    // $ANTLR start "rule__ElkEdgeSection__Group_4_0_0__1"
    // InternalElkGraph.g:3745:1: rule__ElkEdgeSection__Group_4_0_0__1 : rule__ElkEdgeSection__Group_4_0_0__1__Impl rule__ElkEdgeSection__Group_4_0_0__2 ;
    public final void rule__ElkEdgeSection__Group_4_0_0__1() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:3749:1: ( rule__ElkEdgeSection__Group_4_0_0__1__Impl rule__ElkEdgeSection__Group_4_0_0__2 )
            // InternalElkGraph.g:3750:2: rule__ElkEdgeSection__Group_4_0_0__1__Impl rule__ElkEdgeSection__Group_4_0_0__2
            {
            pushFollow(FOLLOW_7);
            rule__ElkEdgeSection__Group_4_0_0__1__Impl();

            state._fsp--;

            pushFollow(FOLLOW_2);
            rule__ElkEdgeSection__Group_4_0_0__2();

            state._fsp--;


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ElkEdgeSection__Group_4_0_0__1"


    // $ANTLR start "rule__ElkEdgeSection__Group_4_0_0__1__Impl"
    // InternalElkGraph.g:3757:1: rule__ElkEdgeSection__Group_4_0_0__1__Impl : ( ':' ) ;
    public final void rule__ElkEdgeSection__Group_4_0_0__1__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:3761:1: ( ( ':' ) )
            // InternalElkGraph.g:3762:1: ( ':' )
            {
            // InternalElkGraph.g:3762:1: ( ':' )
            // InternalElkGraph.g:3763:2: ':'
            {
             before(grammarAccess.getElkEdgeSectionAccess().getColonKeyword_4_0_0_1()); 
            match(input,21,FOLLOW_2); 
             after(grammarAccess.getElkEdgeSectionAccess().getColonKeyword_4_0_0_1()); 

            }


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ElkEdgeSection__Group_4_0_0__1__Impl"


    // $ANTLR start "rule__ElkEdgeSection__Group_4_0_0__2"
    // InternalElkGraph.g:3772:1: rule__ElkEdgeSection__Group_4_0_0__2 : rule__ElkEdgeSection__Group_4_0_0__2__Impl ;
    public final void rule__ElkEdgeSection__Group_4_0_0__2() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:3776:1: ( rule__ElkEdgeSection__Group_4_0_0__2__Impl )
            // InternalElkGraph.g:3777:2: rule__ElkEdgeSection__Group_4_0_0__2__Impl
            {
            pushFollow(FOLLOW_2);
            rule__ElkEdgeSection__Group_4_0_0__2__Impl();

            state._fsp--;


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ElkEdgeSection__Group_4_0_0__2"


    // $ANTLR start "rule__ElkEdgeSection__Group_4_0_0__2__Impl"
    // InternalElkGraph.g:3783:1: rule__ElkEdgeSection__Group_4_0_0__2__Impl : ( ( rule__ElkEdgeSection__IncomingShapeAssignment_4_0_0_2 ) ) ;
    public final void rule__ElkEdgeSection__Group_4_0_0__2__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:3787:1: ( ( ( rule__ElkEdgeSection__IncomingShapeAssignment_4_0_0_2 ) ) )
            // InternalElkGraph.g:3788:1: ( ( rule__ElkEdgeSection__IncomingShapeAssignment_4_0_0_2 ) )
            {
            // InternalElkGraph.g:3788:1: ( ( rule__ElkEdgeSection__IncomingShapeAssignment_4_0_0_2 ) )
            // InternalElkGraph.g:3789:2: ( rule__ElkEdgeSection__IncomingShapeAssignment_4_0_0_2 )
            {
             before(grammarAccess.getElkEdgeSectionAccess().getIncomingShapeAssignment_4_0_0_2()); 
            // InternalElkGraph.g:3790:2: ( rule__ElkEdgeSection__IncomingShapeAssignment_4_0_0_2 )
            // InternalElkGraph.g:3790:3: rule__ElkEdgeSection__IncomingShapeAssignment_4_0_0_2
            {
            pushFollow(FOLLOW_2);
            rule__ElkEdgeSection__IncomingShapeAssignment_4_0_0_2();

            state._fsp--;


            }

             after(grammarAccess.getElkEdgeSectionAccess().getIncomingShapeAssignment_4_0_0_2()); 

            }


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ElkEdgeSection__Group_4_0_0__2__Impl"


    // $ANTLR start "rule__ElkEdgeSection__Group_4_0_1__0"
    // InternalElkGraph.g:3799:1: rule__ElkEdgeSection__Group_4_0_1__0 : rule__ElkEdgeSection__Group_4_0_1__0__Impl rule__ElkEdgeSection__Group_4_0_1__1 ;
    public final void rule__ElkEdgeSection__Group_4_0_1__0() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:3803:1: ( rule__ElkEdgeSection__Group_4_0_1__0__Impl rule__ElkEdgeSection__Group_4_0_1__1 )
            // InternalElkGraph.g:3804:2: rule__ElkEdgeSection__Group_4_0_1__0__Impl rule__ElkEdgeSection__Group_4_0_1__1
            {
            pushFollow(FOLLOW_11);
            rule__ElkEdgeSection__Group_4_0_1__0__Impl();

            state._fsp--;

            pushFollow(FOLLOW_2);
            rule__ElkEdgeSection__Group_4_0_1__1();

            state._fsp--;


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ElkEdgeSection__Group_4_0_1__0"


    // $ANTLR start "rule__ElkEdgeSection__Group_4_0_1__0__Impl"
    // InternalElkGraph.g:3811:1: rule__ElkEdgeSection__Group_4_0_1__0__Impl : ( 'outgoing' ) ;
    public final void rule__ElkEdgeSection__Group_4_0_1__0__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:3815:1: ( ( 'outgoing' ) )
            // InternalElkGraph.g:3816:1: ( 'outgoing' )
            {
            // InternalElkGraph.g:3816:1: ( 'outgoing' )
            // InternalElkGraph.g:3817:2: 'outgoing'
            {
             before(grammarAccess.getElkEdgeSectionAccess().getOutgoingKeyword_4_0_1_0()); 
            match(input,32,FOLLOW_2); 
             after(grammarAccess.getElkEdgeSectionAccess().getOutgoingKeyword_4_0_1_0()); 

            }


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ElkEdgeSection__Group_4_0_1__0__Impl"


    // $ANTLR start "rule__ElkEdgeSection__Group_4_0_1__1"
    // InternalElkGraph.g:3826:1: rule__ElkEdgeSection__Group_4_0_1__1 : rule__ElkEdgeSection__Group_4_0_1__1__Impl rule__ElkEdgeSection__Group_4_0_1__2 ;
    public final void rule__ElkEdgeSection__Group_4_0_1__1() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:3830:1: ( rule__ElkEdgeSection__Group_4_0_1__1__Impl rule__ElkEdgeSection__Group_4_0_1__2 )
            // InternalElkGraph.g:3831:2: rule__ElkEdgeSection__Group_4_0_1__1__Impl rule__ElkEdgeSection__Group_4_0_1__2
            {
            pushFollow(FOLLOW_7);
            rule__ElkEdgeSection__Group_4_0_1__1__Impl();

            state._fsp--;

            pushFollow(FOLLOW_2);
            rule__ElkEdgeSection__Group_4_0_1__2();

            state._fsp--;


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ElkEdgeSection__Group_4_0_1__1"


    // $ANTLR start "rule__ElkEdgeSection__Group_4_0_1__1__Impl"
    // InternalElkGraph.g:3838:1: rule__ElkEdgeSection__Group_4_0_1__1__Impl : ( ':' ) ;
    public final void rule__ElkEdgeSection__Group_4_0_1__1__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:3842:1: ( ( ':' ) )
            // InternalElkGraph.g:3843:1: ( ':' )
            {
            // InternalElkGraph.g:3843:1: ( ':' )
            // InternalElkGraph.g:3844:2: ':'
            {
             before(grammarAccess.getElkEdgeSectionAccess().getColonKeyword_4_0_1_1()); 
            match(input,21,FOLLOW_2); 
             after(grammarAccess.getElkEdgeSectionAccess().getColonKeyword_4_0_1_1()); 

            }


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ElkEdgeSection__Group_4_0_1__1__Impl"


    // $ANTLR start "rule__ElkEdgeSection__Group_4_0_1__2"
    // InternalElkGraph.g:3853:1: rule__ElkEdgeSection__Group_4_0_1__2 : rule__ElkEdgeSection__Group_4_0_1__2__Impl ;
    public final void rule__ElkEdgeSection__Group_4_0_1__2() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:3857:1: ( rule__ElkEdgeSection__Group_4_0_1__2__Impl )
            // InternalElkGraph.g:3858:2: rule__ElkEdgeSection__Group_4_0_1__2__Impl
            {
            pushFollow(FOLLOW_2);
            rule__ElkEdgeSection__Group_4_0_1__2__Impl();

            state._fsp--;


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ElkEdgeSection__Group_4_0_1__2"


    // $ANTLR start "rule__ElkEdgeSection__Group_4_0_1__2__Impl"
    // InternalElkGraph.g:3864:1: rule__ElkEdgeSection__Group_4_0_1__2__Impl : ( ( rule__ElkEdgeSection__OutgoingShapeAssignment_4_0_1_2 ) ) ;
    public final void rule__ElkEdgeSection__Group_4_0_1__2__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:3868:1: ( ( ( rule__ElkEdgeSection__OutgoingShapeAssignment_4_0_1_2 ) ) )
            // InternalElkGraph.g:3869:1: ( ( rule__ElkEdgeSection__OutgoingShapeAssignment_4_0_1_2 ) )
            {
            // InternalElkGraph.g:3869:1: ( ( rule__ElkEdgeSection__OutgoingShapeAssignment_4_0_1_2 ) )
            // InternalElkGraph.g:3870:2: ( rule__ElkEdgeSection__OutgoingShapeAssignment_4_0_1_2 )
            {
             before(grammarAccess.getElkEdgeSectionAccess().getOutgoingShapeAssignment_4_0_1_2()); 
            // InternalElkGraph.g:3871:2: ( rule__ElkEdgeSection__OutgoingShapeAssignment_4_0_1_2 )
            // InternalElkGraph.g:3871:3: rule__ElkEdgeSection__OutgoingShapeAssignment_4_0_1_2
            {
            pushFollow(FOLLOW_2);
            rule__ElkEdgeSection__OutgoingShapeAssignment_4_0_1_2();

            state._fsp--;


            }

             after(grammarAccess.getElkEdgeSectionAccess().getOutgoingShapeAssignment_4_0_1_2()); 

            }


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ElkEdgeSection__Group_4_0_1__2__Impl"


    // $ANTLR start "rule__ElkEdgeSection__Group_4_0_2__0"
    // InternalElkGraph.g:3880:1: rule__ElkEdgeSection__Group_4_0_2__0 : rule__ElkEdgeSection__Group_4_0_2__0__Impl rule__ElkEdgeSection__Group_4_0_2__1 ;
    public final void rule__ElkEdgeSection__Group_4_0_2__0() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:3884:1: ( rule__ElkEdgeSection__Group_4_0_2__0__Impl rule__ElkEdgeSection__Group_4_0_2__1 )
            // InternalElkGraph.g:3885:2: rule__ElkEdgeSection__Group_4_0_2__0__Impl rule__ElkEdgeSection__Group_4_0_2__1
            {
            pushFollow(FOLLOW_11);
            rule__ElkEdgeSection__Group_4_0_2__0__Impl();

            state._fsp--;

            pushFollow(FOLLOW_2);
            rule__ElkEdgeSection__Group_4_0_2__1();

            state._fsp--;


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ElkEdgeSection__Group_4_0_2__0"


    // $ANTLR start "rule__ElkEdgeSection__Group_4_0_2__0__Impl"
    // InternalElkGraph.g:3892:1: rule__ElkEdgeSection__Group_4_0_2__0__Impl : ( 'start' ) ;
    public final void rule__ElkEdgeSection__Group_4_0_2__0__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:3896:1: ( ( 'start' ) )
            // InternalElkGraph.g:3897:1: ( 'start' )
            {
            // InternalElkGraph.g:3897:1: ( 'start' )
            // InternalElkGraph.g:3898:2: 'start'
            {
             before(grammarAccess.getElkEdgeSectionAccess().getStartKeyword_4_0_2_0()); 
            match(input,33,FOLLOW_2); 
             after(grammarAccess.getElkEdgeSectionAccess().getStartKeyword_4_0_2_0()); 

            }


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ElkEdgeSection__Group_4_0_2__0__Impl"


    // $ANTLR start "rule__ElkEdgeSection__Group_4_0_2__1"
    // InternalElkGraph.g:3907:1: rule__ElkEdgeSection__Group_4_0_2__1 : rule__ElkEdgeSection__Group_4_0_2__1__Impl rule__ElkEdgeSection__Group_4_0_2__2 ;
    public final void rule__ElkEdgeSection__Group_4_0_2__1() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:3911:1: ( rule__ElkEdgeSection__Group_4_0_2__1__Impl rule__ElkEdgeSection__Group_4_0_2__2 )
            // InternalElkGraph.g:3912:2: rule__ElkEdgeSection__Group_4_0_2__1__Impl rule__ElkEdgeSection__Group_4_0_2__2
            {
            pushFollow(FOLLOW_17);
            rule__ElkEdgeSection__Group_4_0_2__1__Impl();

            state._fsp--;

            pushFollow(FOLLOW_2);
            rule__ElkEdgeSection__Group_4_0_2__2();

            state._fsp--;


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ElkEdgeSection__Group_4_0_2__1"


    // $ANTLR start "rule__ElkEdgeSection__Group_4_0_2__1__Impl"
    // InternalElkGraph.g:3919:1: rule__ElkEdgeSection__Group_4_0_2__1__Impl : ( ':' ) ;
    public final void rule__ElkEdgeSection__Group_4_0_2__1__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:3923:1: ( ( ':' ) )
            // InternalElkGraph.g:3924:1: ( ':' )
            {
            // InternalElkGraph.g:3924:1: ( ':' )
            // InternalElkGraph.g:3925:2: ':'
            {
             before(grammarAccess.getElkEdgeSectionAccess().getColonKeyword_4_0_2_1()); 
            match(input,21,FOLLOW_2); 
             after(grammarAccess.getElkEdgeSectionAccess().getColonKeyword_4_0_2_1()); 

            }


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ElkEdgeSection__Group_4_0_2__1__Impl"


    // $ANTLR start "rule__ElkEdgeSection__Group_4_0_2__2"
    // InternalElkGraph.g:3934:1: rule__ElkEdgeSection__Group_4_0_2__2 : rule__ElkEdgeSection__Group_4_0_2__2__Impl rule__ElkEdgeSection__Group_4_0_2__3 ;
    public final void rule__ElkEdgeSection__Group_4_0_2__2() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:3938:1: ( rule__ElkEdgeSection__Group_4_0_2__2__Impl rule__ElkEdgeSection__Group_4_0_2__3 )
            // InternalElkGraph.g:3939:2: rule__ElkEdgeSection__Group_4_0_2__2__Impl rule__ElkEdgeSection__Group_4_0_2__3
            {
            pushFollow(FOLLOW_18);
            rule__ElkEdgeSection__Group_4_0_2__2__Impl();

            state._fsp--;

            pushFollow(FOLLOW_2);
            rule__ElkEdgeSection__Group_4_0_2__3();

            state._fsp--;


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ElkEdgeSection__Group_4_0_2__2"


    // $ANTLR start "rule__ElkEdgeSection__Group_4_0_2__2__Impl"
    // InternalElkGraph.g:3946:1: rule__ElkEdgeSection__Group_4_0_2__2__Impl : ( ( rule__ElkEdgeSection__StartXAssignment_4_0_2_2 ) ) ;
    public final void rule__ElkEdgeSection__Group_4_0_2__2__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:3950:1: ( ( ( rule__ElkEdgeSection__StartXAssignment_4_0_2_2 ) ) )
            // InternalElkGraph.g:3951:1: ( ( rule__ElkEdgeSection__StartXAssignment_4_0_2_2 ) )
            {
            // InternalElkGraph.g:3951:1: ( ( rule__ElkEdgeSection__StartXAssignment_4_0_2_2 ) )
            // InternalElkGraph.g:3952:2: ( rule__ElkEdgeSection__StartXAssignment_4_0_2_2 )
            {
             before(grammarAccess.getElkEdgeSectionAccess().getStartXAssignment_4_0_2_2()); 
            // InternalElkGraph.g:3953:2: ( rule__ElkEdgeSection__StartXAssignment_4_0_2_2 )
            // InternalElkGraph.g:3953:3: rule__ElkEdgeSection__StartXAssignment_4_0_2_2
            {
            pushFollow(FOLLOW_2);
            rule__ElkEdgeSection__StartXAssignment_4_0_2_2();

            state._fsp--;


            }

             after(grammarAccess.getElkEdgeSectionAccess().getStartXAssignment_4_0_2_2()); 

            }


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ElkEdgeSection__Group_4_0_2__2__Impl"


    // $ANTLR start "rule__ElkEdgeSection__Group_4_0_2__3"
    // InternalElkGraph.g:3961:1: rule__ElkEdgeSection__Group_4_0_2__3 : rule__ElkEdgeSection__Group_4_0_2__3__Impl rule__ElkEdgeSection__Group_4_0_2__4 ;
    public final void rule__ElkEdgeSection__Group_4_0_2__3() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:3965:1: ( rule__ElkEdgeSection__Group_4_0_2__3__Impl rule__ElkEdgeSection__Group_4_0_2__4 )
            // InternalElkGraph.g:3966:2: rule__ElkEdgeSection__Group_4_0_2__3__Impl rule__ElkEdgeSection__Group_4_0_2__4
            {
            pushFollow(FOLLOW_17);
            rule__ElkEdgeSection__Group_4_0_2__3__Impl();

            state._fsp--;

            pushFollow(FOLLOW_2);
            rule__ElkEdgeSection__Group_4_0_2__4();

            state._fsp--;


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ElkEdgeSection__Group_4_0_2__3"


    // $ANTLR start "rule__ElkEdgeSection__Group_4_0_2__3__Impl"
    // InternalElkGraph.g:3973:1: rule__ElkEdgeSection__Group_4_0_2__3__Impl : ( ',' ) ;
    public final void rule__ElkEdgeSection__Group_4_0_2__3__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:3977:1: ( ( ',' ) )
            // InternalElkGraph.g:3978:1: ( ',' )
            {
            // InternalElkGraph.g:3978:1: ( ',' )
            // InternalElkGraph.g:3979:2: ','
            {
             before(grammarAccess.getElkEdgeSectionAccess().getCommaKeyword_4_0_2_3()); 
            match(input,27,FOLLOW_2); 
             after(grammarAccess.getElkEdgeSectionAccess().getCommaKeyword_4_0_2_3()); 

            }


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ElkEdgeSection__Group_4_0_2__3__Impl"


    // $ANTLR start "rule__ElkEdgeSection__Group_4_0_2__4"
    // InternalElkGraph.g:3988:1: rule__ElkEdgeSection__Group_4_0_2__4 : rule__ElkEdgeSection__Group_4_0_2__4__Impl ;
    public final void rule__ElkEdgeSection__Group_4_0_2__4() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:3992:1: ( rule__ElkEdgeSection__Group_4_0_2__4__Impl )
            // InternalElkGraph.g:3993:2: rule__ElkEdgeSection__Group_4_0_2__4__Impl
            {
            pushFollow(FOLLOW_2);
            rule__ElkEdgeSection__Group_4_0_2__4__Impl();

            state._fsp--;


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ElkEdgeSection__Group_4_0_2__4"


    // $ANTLR start "rule__ElkEdgeSection__Group_4_0_2__4__Impl"
    // InternalElkGraph.g:3999:1: rule__ElkEdgeSection__Group_4_0_2__4__Impl : ( ( rule__ElkEdgeSection__StartYAssignment_4_0_2_4 ) ) ;
    public final void rule__ElkEdgeSection__Group_4_0_2__4__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:4003:1: ( ( ( rule__ElkEdgeSection__StartYAssignment_4_0_2_4 ) ) )
            // InternalElkGraph.g:4004:1: ( ( rule__ElkEdgeSection__StartYAssignment_4_0_2_4 ) )
            {
            // InternalElkGraph.g:4004:1: ( ( rule__ElkEdgeSection__StartYAssignment_4_0_2_4 ) )
            // InternalElkGraph.g:4005:2: ( rule__ElkEdgeSection__StartYAssignment_4_0_2_4 )
            {
             before(grammarAccess.getElkEdgeSectionAccess().getStartYAssignment_4_0_2_4()); 
            // InternalElkGraph.g:4006:2: ( rule__ElkEdgeSection__StartYAssignment_4_0_2_4 )
            // InternalElkGraph.g:4006:3: rule__ElkEdgeSection__StartYAssignment_4_0_2_4
            {
            pushFollow(FOLLOW_2);
            rule__ElkEdgeSection__StartYAssignment_4_0_2_4();

            state._fsp--;


            }

             after(grammarAccess.getElkEdgeSectionAccess().getStartYAssignment_4_0_2_4()); 

            }


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ElkEdgeSection__Group_4_0_2__4__Impl"


    // $ANTLR start "rule__ElkEdgeSection__Group_4_0_3__0"
    // InternalElkGraph.g:4015:1: rule__ElkEdgeSection__Group_4_0_3__0 : rule__ElkEdgeSection__Group_4_0_3__0__Impl rule__ElkEdgeSection__Group_4_0_3__1 ;
    public final void rule__ElkEdgeSection__Group_4_0_3__0() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:4019:1: ( rule__ElkEdgeSection__Group_4_0_3__0__Impl rule__ElkEdgeSection__Group_4_0_3__1 )
            // InternalElkGraph.g:4020:2: rule__ElkEdgeSection__Group_4_0_3__0__Impl rule__ElkEdgeSection__Group_4_0_3__1
            {
            pushFollow(FOLLOW_11);
            rule__ElkEdgeSection__Group_4_0_3__0__Impl();

            state._fsp--;

            pushFollow(FOLLOW_2);
            rule__ElkEdgeSection__Group_4_0_3__1();

            state._fsp--;


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ElkEdgeSection__Group_4_0_3__0"


    // $ANTLR start "rule__ElkEdgeSection__Group_4_0_3__0__Impl"
    // InternalElkGraph.g:4027:1: rule__ElkEdgeSection__Group_4_0_3__0__Impl : ( 'end' ) ;
    public final void rule__ElkEdgeSection__Group_4_0_3__0__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:4031:1: ( ( 'end' ) )
            // InternalElkGraph.g:4032:1: ( 'end' )
            {
            // InternalElkGraph.g:4032:1: ( 'end' )
            // InternalElkGraph.g:4033:2: 'end'
            {
             before(grammarAccess.getElkEdgeSectionAccess().getEndKeyword_4_0_3_0()); 
            match(input,34,FOLLOW_2); 
             after(grammarAccess.getElkEdgeSectionAccess().getEndKeyword_4_0_3_0()); 

            }


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ElkEdgeSection__Group_4_0_3__0__Impl"


    // $ANTLR start "rule__ElkEdgeSection__Group_4_0_3__1"
    // InternalElkGraph.g:4042:1: rule__ElkEdgeSection__Group_4_0_3__1 : rule__ElkEdgeSection__Group_4_0_3__1__Impl rule__ElkEdgeSection__Group_4_0_3__2 ;
    public final void rule__ElkEdgeSection__Group_4_0_3__1() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:4046:1: ( rule__ElkEdgeSection__Group_4_0_3__1__Impl rule__ElkEdgeSection__Group_4_0_3__2 )
            // InternalElkGraph.g:4047:2: rule__ElkEdgeSection__Group_4_0_3__1__Impl rule__ElkEdgeSection__Group_4_0_3__2
            {
            pushFollow(FOLLOW_17);
            rule__ElkEdgeSection__Group_4_0_3__1__Impl();

            state._fsp--;

            pushFollow(FOLLOW_2);
            rule__ElkEdgeSection__Group_4_0_3__2();

            state._fsp--;


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ElkEdgeSection__Group_4_0_3__1"


    // $ANTLR start "rule__ElkEdgeSection__Group_4_0_3__1__Impl"
    // InternalElkGraph.g:4054:1: rule__ElkEdgeSection__Group_4_0_3__1__Impl : ( ':' ) ;
    public final void rule__ElkEdgeSection__Group_4_0_3__1__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:4058:1: ( ( ':' ) )
            // InternalElkGraph.g:4059:1: ( ':' )
            {
            // InternalElkGraph.g:4059:1: ( ':' )
            // InternalElkGraph.g:4060:2: ':'
            {
             before(grammarAccess.getElkEdgeSectionAccess().getColonKeyword_4_0_3_1()); 
            match(input,21,FOLLOW_2); 
             after(grammarAccess.getElkEdgeSectionAccess().getColonKeyword_4_0_3_1()); 

            }


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ElkEdgeSection__Group_4_0_3__1__Impl"


    // $ANTLR start "rule__ElkEdgeSection__Group_4_0_3__2"
    // InternalElkGraph.g:4069:1: rule__ElkEdgeSection__Group_4_0_3__2 : rule__ElkEdgeSection__Group_4_0_3__2__Impl rule__ElkEdgeSection__Group_4_0_3__3 ;
    public final void rule__ElkEdgeSection__Group_4_0_3__2() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:4073:1: ( rule__ElkEdgeSection__Group_4_0_3__2__Impl rule__ElkEdgeSection__Group_4_0_3__3 )
            // InternalElkGraph.g:4074:2: rule__ElkEdgeSection__Group_4_0_3__2__Impl rule__ElkEdgeSection__Group_4_0_3__3
            {
            pushFollow(FOLLOW_18);
            rule__ElkEdgeSection__Group_4_0_3__2__Impl();

            state._fsp--;

            pushFollow(FOLLOW_2);
            rule__ElkEdgeSection__Group_4_0_3__3();

            state._fsp--;


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ElkEdgeSection__Group_4_0_3__2"


    // $ANTLR start "rule__ElkEdgeSection__Group_4_0_3__2__Impl"
    // InternalElkGraph.g:4081:1: rule__ElkEdgeSection__Group_4_0_3__2__Impl : ( ( rule__ElkEdgeSection__EndXAssignment_4_0_3_2 ) ) ;
    public final void rule__ElkEdgeSection__Group_4_0_3__2__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:4085:1: ( ( ( rule__ElkEdgeSection__EndXAssignment_4_0_3_2 ) ) )
            // InternalElkGraph.g:4086:1: ( ( rule__ElkEdgeSection__EndXAssignment_4_0_3_2 ) )
            {
            // InternalElkGraph.g:4086:1: ( ( rule__ElkEdgeSection__EndXAssignment_4_0_3_2 ) )
            // InternalElkGraph.g:4087:2: ( rule__ElkEdgeSection__EndXAssignment_4_0_3_2 )
            {
             before(grammarAccess.getElkEdgeSectionAccess().getEndXAssignment_4_0_3_2()); 
            // InternalElkGraph.g:4088:2: ( rule__ElkEdgeSection__EndXAssignment_4_0_3_2 )
            // InternalElkGraph.g:4088:3: rule__ElkEdgeSection__EndXAssignment_4_0_3_2
            {
            pushFollow(FOLLOW_2);
            rule__ElkEdgeSection__EndXAssignment_4_0_3_2();

            state._fsp--;


            }

             after(grammarAccess.getElkEdgeSectionAccess().getEndXAssignment_4_0_3_2()); 

            }


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ElkEdgeSection__Group_4_0_3__2__Impl"


    // $ANTLR start "rule__ElkEdgeSection__Group_4_0_3__3"
    // InternalElkGraph.g:4096:1: rule__ElkEdgeSection__Group_4_0_3__3 : rule__ElkEdgeSection__Group_4_0_3__3__Impl rule__ElkEdgeSection__Group_4_0_3__4 ;
    public final void rule__ElkEdgeSection__Group_4_0_3__3() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:4100:1: ( rule__ElkEdgeSection__Group_4_0_3__3__Impl rule__ElkEdgeSection__Group_4_0_3__4 )
            // InternalElkGraph.g:4101:2: rule__ElkEdgeSection__Group_4_0_3__3__Impl rule__ElkEdgeSection__Group_4_0_3__4
            {
            pushFollow(FOLLOW_17);
            rule__ElkEdgeSection__Group_4_0_3__3__Impl();

            state._fsp--;

            pushFollow(FOLLOW_2);
            rule__ElkEdgeSection__Group_4_0_3__4();

            state._fsp--;


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ElkEdgeSection__Group_4_0_3__3"


    // $ANTLR start "rule__ElkEdgeSection__Group_4_0_3__3__Impl"
    // InternalElkGraph.g:4108:1: rule__ElkEdgeSection__Group_4_0_3__3__Impl : ( ',' ) ;
    public final void rule__ElkEdgeSection__Group_4_0_3__3__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:4112:1: ( ( ',' ) )
            // InternalElkGraph.g:4113:1: ( ',' )
            {
            // InternalElkGraph.g:4113:1: ( ',' )
            // InternalElkGraph.g:4114:2: ','
            {
             before(grammarAccess.getElkEdgeSectionAccess().getCommaKeyword_4_0_3_3()); 
            match(input,27,FOLLOW_2); 
             after(grammarAccess.getElkEdgeSectionAccess().getCommaKeyword_4_0_3_3()); 

            }


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ElkEdgeSection__Group_4_0_3__3__Impl"


    // $ANTLR start "rule__ElkEdgeSection__Group_4_0_3__4"
    // InternalElkGraph.g:4123:1: rule__ElkEdgeSection__Group_4_0_3__4 : rule__ElkEdgeSection__Group_4_0_3__4__Impl ;
    public final void rule__ElkEdgeSection__Group_4_0_3__4() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:4127:1: ( rule__ElkEdgeSection__Group_4_0_3__4__Impl )
            // InternalElkGraph.g:4128:2: rule__ElkEdgeSection__Group_4_0_3__4__Impl
            {
            pushFollow(FOLLOW_2);
            rule__ElkEdgeSection__Group_4_0_3__4__Impl();

            state._fsp--;


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ElkEdgeSection__Group_4_0_3__4"


    // $ANTLR start "rule__ElkEdgeSection__Group_4_0_3__4__Impl"
    // InternalElkGraph.g:4134:1: rule__ElkEdgeSection__Group_4_0_3__4__Impl : ( ( rule__ElkEdgeSection__EndYAssignment_4_0_3_4 ) ) ;
    public final void rule__ElkEdgeSection__Group_4_0_3__4__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:4138:1: ( ( ( rule__ElkEdgeSection__EndYAssignment_4_0_3_4 ) ) )
            // InternalElkGraph.g:4139:1: ( ( rule__ElkEdgeSection__EndYAssignment_4_0_3_4 ) )
            {
            // InternalElkGraph.g:4139:1: ( ( rule__ElkEdgeSection__EndYAssignment_4_0_3_4 ) )
            // InternalElkGraph.g:4140:2: ( rule__ElkEdgeSection__EndYAssignment_4_0_3_4 )
            {
             before(grammarAccess.getElkEdgeSectionAccess().getEndYAssignment_4_0_3_4()); 
            // InternalElkGraph.g:4141:2: ( rule__ElkEdgeSection__EndYAssignment_4_0_3_4 )
            // InternalElkGraph.g:4141:3: rule__ElkEdgeSection__EndYAssignment_4_0_3_4
            {
            pushFollow(FOLLOW_2);
            rule__ElkEdgeSection__EndYAssignment_4_0_3_4();

            state._fsp--;


            }

             after(grammarAccess.getElkEdgeSectionAccess().getEndYAssignment_4_0_3_4()); 

            }


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ElkEdgeSection__Group_4_0_3__4__Impl"


    // $ANTLR start "rule__ElkEdgeSection__Group_4_1__0"
    // InternalElkGraph.g:4150:1: rule__ElkEdgeSection__Group_4_1__0 : rule__ElkEdgeSection__Group_4_1__0__Impl rule__ElkEdgeSection__Group_4_1__1 ;
    public final void rule__ElkEdgeSection__Group_4_1__0() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:4154:1: ( rule__ElkEdgeSection__Group_4_1__0__Impl rule__ElkEdgeSection__Group_4_1__1 )
            // InternalElkGraph.g:4155:2: rule__ElkEdgeSection__Group_4_1__0__Impl rule__ElkEdgeSection__Group_4_1__1
            {
            pushFollow(FOLLOW_11);
            rule__ElkEdgeSection__Group_4_1__0__Impl();

            state._fsp--;

            pushFollow(FOLLOW_2);
            rule__ElkEdgeSection__Group_4_1__1();

            state._fsp--;


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ElkEdgeSection__Group_4_1__0"


    // $ANTLR start "rule__ElkEdgeSection__Group_4_1__0__Impl"
    // InternalElkGraph.g:4162:1: rule__ElkEdgeSection__Group_4_1__0__Impl : ( 'bends' ) ;
    public final void rule__ElkEdgeSection__Group_4_1__0__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:4166:1: ( ( 'bends' ) )
            // InternalElkGraph.g:4167:1: ( 'bends' )
            {
            // InternalElkGraph.g:4167:1: ( 'bends' )
            // InternalElkGraph.g:4168:2: 'bends'
            {
             before(grammarAccess.getElkEdgeSectionAccess().getBendsKeyword_4_1_0()); 
            match(input,35,FOLLOW_2); 
             after(grammarAccess.getElkEdgeSectionAccess().getBendsKeyword_4_1_0()); 

            }


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ElkEdgeSection__Group_4_1__0__Impl"


    // $ANTLR start "rule__ElkEdgeSection__Group_4_1__1"
    // InternalElkGraph.g:4177:1: rule__ElkEdgeSection__Group_4_1__1 : rule__ElkEdgeSection__Group_4_1__1__Impl rule__ElkEdgeSection__Group_4_1__2 ;
    public final void rule__ElkEdgeSection__Group_4_1__1() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:4181:1: ( rule__ElkEdgeSection__Group_4_1__1__Impl rule__ElkEdgeSection__Group_4_1__2 )
            // InternalElkGraph.g:4182:2: rule__ElkEdgeSection__Group_4_1__1__Impl rule__ElkEdgeSection__Group_4_1__2
            {
            pushFollow(FOLLOW_17);
            rule__ElkEdgeSection__Group_4_1__1__Impl();

            state._fsp--;

            pushFollow(FOLLOW_2);
            rule__ElkEdgeSection__Group_4_1__2();

            state._fsp--;


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ElkEdgeSection__Group_4_1__1"


    // $ANTLR start "rule__ElkEdgeSection__Group_4_1__1__Impl"
    // InternalElkGraph.g:4189:1: rule__ElkEdgeSection__Group_4_1__1__Impl : ( ':' ) ;
    public final void rule__ElkEdgeSection__Group_4_1__1__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:4193:1: ( ( ':' ) )
            // InternalElkGraph.g:4194:1: ( ':' )
            {
            // InternalElkGraph.g:4194:1: ( ':' )
            // InternalElkGraph.g:4195:2: ':'
            {
             before(grammarAccess.getElkEdgeSectionAccess().getColonKeyword_4_1_1()); 
            match(input,21,FOLLOW_2); 
             after(grammarAccess.getElkEdgeSectionAccess().getColonKeyword_4_1_1()); 

            }


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ElkEdgeSection__Group_4_1__1__Impl"


    // $ANTLR start "rule__ElkEdgeSection__Group_4_1__2"
    // InternalElkGraph.g:4204:1: rule__ElkEdgeSection__Group_4_1__2 : rule__ElkEdgeSection__Group_4_1__2__Impl rule__ElkEdgeSection__Group_4_1__3 ;
    public final void rule__ElkEdgeSection__Group_4_1__2() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:4208:1: ( rule__ElkEdgeSection__Group_4_1__2__Impl rule__ElkEdgeSection__Group_4_1__3 )
            // InternalElkGraph.g:4209:2: rule__ElkEdgeSection__Group_4_1__2__Impl rule__ElkEdgeSection__Group_4_1__3
            {
            pushFollow(FOLLOW_25);
            rule__ElkEdgeSection__Group_4_1__2__Impl();

            state._fsp--;

            pushFollow(FOLLOW_2);
            rule__ElkEdgeSection__Group_4_1__3();

            state._fsp--;


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ElkEdgeSection__Group_4_1__2"


    // $ANTLR start "rule__ElkEdgeSection__Group_4_1__2__Impl"
    // InternalElkGraph.g:4216:1: rule__ElkEdgeSection__Group_4_1__2__Impl : ( ( rule__ElkEdgeSection__BendPointsAssignment_4_1_2 ) ) ;
    public final void rule__ElkEdgeSection__Group_4_1__2__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:4220:1: ( ( ( rule__ElkEdgeSection__BendPointsAssignment_4_1_2 ) ) )
            // InternalElkGraph.g:4221:1: ( ( rule__ElkEdgeSection__BendPointsAssignment_4_1_2 ) )
            {
            // InternalElkGraph.g:4221:1: ( ( rule__ElkEdgeSection__BendPointsAssignment_4_1_2 ) )
            // InternalElkGraph.g:4222:2: ( rule__ElkEdgeSection__BendPointsAssignment_4_1_2 )
            {
             before(grammarAccess.getElkEdgeSectionAccess().getBendPointsAssignment_4_1_2()); 
            // InternalElkGraph.g:4223:2: ( rule__ElkEdgeSection__BendPointsAssignment_4_1_2 )
            // InternalElkGraph.g:4223:3: rule__ElkEdgeSection__BendPointsAssignment_4_1_2
            {
            pushFollow(FOLLOW_2);
            rule__ElkEdgeSection__BendPointsAssignment_4_1_2();

            state._fsp--;


            }

             after(grammarAccess.getElkEdgeSectionAccess().getBendPointsAssignment_4_1_2()); 

            }


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ElkEdgeSection__Group_4_1__2__Impl"


    // $ANTLR start "rule__ElkEdgeSection__Group_4_1__3"
    // InternalElkGraph.g:4231:1: rule__ElkEdgeSection__Group_4_1__3 : rule__ElkEdgeSection__Group_4_1__3__Impl ;
    public final void rule__ElkEdgeSection__Group_4_1__3() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:4235:1: ( rule__ElkEdgeSection__Group_4_1__3__Impl )
            // InternalElkGraph.g:4236:2: rule__ElkEdgeSection__Group_4_1__3__Impl
            {
            pushFollow(FOLLOW_2);
            rule__ElkEdgeSection__Group_4_1__3__Impl();

            state._fsp--;


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ElkEdgeSection__Group_4_1__3"


    // $ANTLR start "rule__ElkEdgeSection__Group_4_1__3__Impl"
    // InternalElkGraph.g:4242:1: rule__ElkEdgeSection__Group_4_1__3__Impl : ( ( rule__ElkEdgeSection__Group_4_1_3__0 )* ) ;
    public final void rule__ElkEdgeSection__Group_4_1__3__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:4246:1: ( ( ( rule__ElkEdgeSection__Group_4_1_3__0 )* ) )
            // InternalElkGraph.g:4247:1: ( ( rule__ElkEdgeSection__Group_4_1_3__0 )* )
            {
            // InternalElkGraph.g:4247:1: ( ( rule__ElkEdgeSection__Group_4_1_3__0 )* )
            // InternalElkGraph.g:4248:2: ( rule__ElkEdgeSection__Group_4_1_3__0 )*
            {
             before(grammarAccess.getElkEdgeSectionAccess().getGroup_4_1_3()); 
            // InternalElkGraph.g:4249:2: ( rule__ElkEdgeSection__Group_4_1_3__0 )*
            loop40:
            do {
                int alt40=2;
                int LA40_0 = input.LA(1);

                if ( (LA40_0==36) ) {
                    alt40=1;
                }


                switch (alt40) {
            	case 1 :
            	    // InternalElkGraph.g:4249:3: rule__ElkEdgeSection__Group_4_1_3__0
            	    {
            	    pushFollow(FOLLOW_26);
            	    rule__ElkEdgeSection__Group_4_1_3__0();

            	    state._fsp--;


            	    }
            	    break;

            	default :
            	    break loop40;
                }
            } while (true);

             after(grammarAccess.getElkEdgeSectionAccess().getGroup_4_1_3()); 

            }


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ElkEdgeSection__Group_4_1__3__Impl"


    // $ANTLR start "rule__ElkEdgeSection__Group_4_1_3__0"
    // InternalElkGraph.g:4258:1: rule__ElkEdgeSection__Group_4_1_3__0 : rule__ElkEdgeSection__Group_4_1_3__0__Impl rule__ElkEdgeSection__Group_4_1_3__1 ;
    public final void rule__ElkEdgeSection__Group_4_1_3__0() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:4262:1: ( rule__ElkEdgeSection__Group_4_1_3__0__Impl rule__ElkEdgeSection__Group_4_1_3__1 )
            // InternalElkGraph.g:4263:2: rule__ElkEdgeSection__Group_4_1_3__0__Impl rule__ElkEdgeSection__Group_4_1_3__1
            {
            pushFollow(FOLLOW_17);
            rule__ElkEdgeSection__Group_4_1_3__0__Impl();

            state._fsp--;

            pushFollow(FOLLOW_2);
            rule__ElkEdgeSection__Group_4_1_3__1();

            state._fsp--;


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ElkEdgeSection__Group_4_1_3__0"


    // $ANTLR start "rule__ElkEdgeSection__Group_4_1_3__0__Impl"
    // InternalElkGraph.g:4270:1: rule__ElkEdgeSection__Group_4_1_3__0__Impl : ( '|' ) ;
    public final void rule__ElkEdgeSection__Group_4_1_3__0__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:4274:1: ( ( '|' ) )
            // InternalElkGraph.g:4275:1: ( '|' )
            {
            // InternalElkGraph.g:4275:1: ( '|' )
            // InternalElkGraph.g:4276:2: '|'
            {
             before(grammarAccess.getElkEdgeSectionAccess().getVerticalLineKeyword_4_1_3_0()); 
            match(input,36,FOLLOW_2); 
             after(grammarAccess.getElkEdgeSectionAccess().getVerticalLineKeyword_4_1_3_0()); 

            }


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ElkEdgeSection__Group_4_1_3__0__Impl"


    // $ANTLR start "rule__ElkEdgeSection__Group_4_1_3__1"
    // InternalElkGraph.g:4285:1: rule__ElkEdgeSection__Group_4_1_3__1 : rule__ElkEdgeSection__Group_4_1_3__1__Impl ;
    public final void rule__ElkEdgeSection__Group_4_1_3__1() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:4289:1: ( rule__ElkEdgeSection__Group_4_1_3__1__Impl )
            // InternalElkGraph.g:4290:2: rule__ElkEdgeSection__Group_4_1_3__1__Impl
            {
            pushFollow(FOLLOW_2);
            rule__ElkEdgeSection__Group_4_1_3__1__Impl();

            state._fsp--;


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ElkEdgeSection__Group_4_1_3__1"


    // $ANTLR start "rule__ElkEdgeSection__Group_4_1_3__1__Impl"
    // InternalElkGraph.g:4296:1: rule__ElkEdgeSection__Group_4_1_3__1__Impl : ( ( rule__ElkEdgeSection__BendPointsAssignment_4_1_3_1 ) ) ;
    public final void rule__ElkEdgeSection__Group_4_1_3__1__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:4300:1: ( ( ( rule__ElkEdgeSection__BendPointsAssignment_4_1_3_1 ) ) )
            // InternalElkGraph.g:4301:1: ( ( rule__ElkEdgeSection__BendPointsAssignment_4_1_3_1 ) )
            {
            // InternalElkGraph.g:4301:1: ( ( rule__ElkEdgeSection__BendPointsAssignment_4_1_3_1 ) )
            // InternalElkGraph.g:4302:2: ( rule__ElkEdgeSection__BendPointsAssignment_4_1_3_1 )
            {
             before(grammarAccess.getElkEdgeSectionAccess().getBendPointsAssignment_4_1_3_1()); 
            // InternalElkGraph.g:4303:2: ( rule__ElkEdgeSection__BendPointsAssignment_4_1_3_1 )
            // InternalElkGraph.g:4303:3: rule__ElkEdgeSection__BendPointsAssignment_4_1_3_1
            {
            pushFollow(FOLLOW_2);
            rule__ElkEdgeSection__BendPointsAssignment_4_1_3_1();

            state._fsp--;


            }

             after(grammarAccess.getElkEdgeSectionAccess().getBendPointsAssignment_4_1_3_1()); 

            }


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ElkEdgeSection__Group_4_1_3__1__Impl"


    // $ANTLR start "rule__ElkBendPoint__Group__0"
    // InternalElkGraph.g:4312:1: rule__ElkBendPoint__Group__0 : rule__ElkBendPoint__Group__0__Impl rule__ElkBendPoint__Group__1 ;
    public final void rule__ElkBendPoint__Group__0() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:4316:1: ( rule__ElkBendPoint__Group__0__Impl rule__ElkBendPoint__Group__1 )
            // InternalElkGraph.g:4317:2: rule__ElkBendPoint__Group__0__Impl rule__ElkBendPoint__Group__1
            {
            pushFollow(FOLLOW_18);
            rule__ElkBendPoint__Group__0__Impl();

            state._fsp--;

            pushFollow(FOLLOW_2);
            rule__ElkBendPoint__Group__1();

            state._fsp--;


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ElkBendPoint__Group__0"


    // $ANTLR start "rule__ElkBendPoint__Group__0__Impl"
    // InternalElkGraph.g:4324:1: rule__ElkBendPoint__Group__0__Impl : ( ( rule__ElkBendPoint__XAssignment_0 ) ) ;
    public final void rule__ElkBendPoint__Group__0__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:4328:1: ( ( ( rule__ElkBendPoint__XAssignment_0 ) ) )
            // InternalElkGraph.g:4329:1: ( ( rule__ElkBendPoint__XAssignment_0 ) )
            {
            // InternalElkGraph.g:4329:1: ( ( rule__ElkBendPoint__XAssignment_0 ) )
            // InternalElkGraph.g:4330:2: ( rule__ElkBendPoint__XAssignment_0 )
            {
             before(grammarAccess.getElkBendPointAccess().getXAssignment_0()); 
            // InternalElkGraph.g:4331:2: ( rule__ElkBendPoint__XAssignment_0 )
            // InternalElkGraph.g:4331:3: rule__ElkBendPoint__XAssignment_0
            {
            pushFollow(FOLLOW_2);
            rule__ElkBendPoint__XAssignment_0();

            state._fsp--;


            }

             after(grammarAccess.getElkBendPointAccess().getXAssignment_0()); 

            }


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ElkBendPoint__Group__0__Impl"


    // $ANTLR start "rule__ElkBendPoint__Group__1"
    // InternalElkGraph.g:4339:1: rule__ElkBendPoint__Group__1 : rule__ElkBendPoint__Group__1__Impl rule__ElkBendPoint__Group__2 ;
    public final void rule__ElkBendPoint__Group__1() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:4343:1: ( rule__ElkBendPoint__Group__1__Impl rule__ElkBendPoint__Group__2 )
            // InternalElkGraph.g:4344:2: rule__ElkBendPoint__Group__1__Impl rule__ElkBendPoint__Group__2
            {
            pushFollow(FOLLOW_17);
            rule__ElkBendPoint__Group__1__Impl();

            state._fsp--;

            pushFollow(FOLLOW_2);
            rule__ElkBendPoint__Group__2();

            state._fsp--;


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ElkBendPoint__Group__1"


    // $ANTLR start "rule__ElkBendPoint__Group__1__Impl"
    // InternalElkGraph.g:4351:1: rule__ElkBendPoint__Group__1__Impl : ( ',' ) ;
    public final void rule__ElkBendPoint__Group__1__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:4355:1: ( ( ',' ) )
            // InternalElkGraph.g:4356:1: ( ',' )
            {
            // InternalElkGraph.g:4356:1: ( ',' )
            // InternalElkGraph.g:4357:2: ','
            {
             before(grammarAccess.getElkBendPointAccess().getCommaKeyword_1()); 
            match(input,27,FOLLOW_2); 
             after(grammarAccess.getElkBendPointAccess().getCommaKeyword_1()); 

            }


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ElkBendPoint__Group__1__Impl"


    // $ANTLR start "rule__ElkBendPoint__Group__2"
    // InternalElkGraph.g:4366:1: rule__ElkBendPoint__Group__2 : rule__ElkBendPoint__Group__2__Impl ;
    public final void rule__ElkBendPoint__Group__2() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:4370:1: ( rule__ElkBendPoint__Group__2__Impl )
            // InternalElkGraph.g:4371:2: rule__ElkBendPoint__Group__2__Impl
            {
            pushFollow(FOLLOW_2);
            rule__ElkBendPoint__Group__2__Impl();

            state._fsp--;


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ElkBendPoint__Group__2"


    // $ANTLR start "rule__ElkBendPoint__Group__2__Impl"
    // InternalElkGraph.g:4377:1: rule__ElkBendPoint__Group__2__Impl : ( ( rule__ElkBendPoint__YAssignment_2 ) ) ;
    public final void rule__ElkBendPoint__Group__2__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:4381:1: ( ( ( rule__ElkBendPoint__YAssignment_2 ) ) )
            // InternalElkGraph.g:4382:1: ( ( rule__ElkBendPoint__YAssignment_2 ) )
            {
            // InternalElkGraph.g:4382:1: ( ( rule__ElkBendPoint__YAssignment_2 ) )
            // InternalElkGraph.g:4383:2: ( rule__ElkBendPoint__YAssignment_2 )
            {
             before(grammarAccess.getElkBendPointAccess().getYAssignment_2()); 
            // InternalElkGraph.g:4384:2: ( rule__ElkBendPoint__YAssignment_2 )
            // InternalElkGraph.g:4384:3: rule__ElkBendPoint__YAssignment_2
            {
            pushFollow(FOLLOW_2);
            rule__ElkBendPoint__YAssignment_2();

            state._fsp--;


            }

             after(grammarAccess.getElkBendPointAccess().getYAssignment_2()); 

            }


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ElkBendPoint__Group__2__Impl"


    // $ANTLR start "rule__QualifiedId__Group__0"
    // InternalElkGraph.g:4393:1: rule__QualifiedId__Group__0 : rule__QualifiedId__Group__0__Impl rule__QualifiedId__Group__1 ;
    public final void rule__QualifiedId__Group__0() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:4397:1: ( rule__QualifiedId__Group__0__Impl rule__QualifiedId__Group__1 )
            // InternalElkGraph.g:4398:2: rule__QualifiedId__Group__0__Impl rule__QualifiedId__Group__1
            {
            pushFollow(FOLLOW_28);
            rule__QualifiedId__Group__0__Impl();

            state._fsp--;

            pushFollow(FOLLOW_2);
            rule__QualifiedId__Group__1();

            state._fsp--;


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__QualifiedId__Group__0"


    // $ANTLR start "rule__QualifiedId__Group__0__Impl"
    // InternalElkGraph.g:4405:1: rule__QualifiedId__Group__0__Impl : ( RULE_ID ) ;
    public final void rule__QualifiedId__Group__0__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:4409:1: ( ( RULE_ID ) )
            // InternalElkGraph.g:4410:1: ( RULE_ID )
            {
            // InternalElkGraph.g:4410:1: ( RULE_ID )
            // InternalElkGraph.g:4411:2: RULE_ID
            {
             before(grammarAccess.getQualifiedIdAccess().getIDTerminalRuleCall_0()); 
            match(input,RULE_ID,FOLLOW_2); 
             after(grammarAccess.getQualifiedIdAccess().getIDTerminalRuleCall_0()); 

            }


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__QualifiedId__Group__0__Impl"


    // $ANTLR start "rule__QualifiedId__Group__1"
    // InternalElkGraph.g:4420:1: rule__QualifiedId__Group__1 : rule__QualifiedId__Group__1__Impl ;
    public final void rule__QualifiedId__Group__1() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:4424:1: ( rule__QualifiedId__Group__1__Impl )
            // InternalElkGraph.g:4425:2: rule__QualifiedId__Group__1__Impl
            {
            pushFollow(FOLLOW_2);
            rule__QualifiedId__Group__1__Impl();

            state._fsp--;


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__QualifiedId__Group__1"


    // $ANTLR start "rule__QualifiedId__Group__1__Impl"
    // InternalElkGraph.g:4431:1: rule__QualifiedId__Group__1__Impl : ( ( rule__QualifiedId__Group_1__0 )* ) ;
    public final void rule__QualifiedId__Group__1__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:4435:1: ( ( ( rule__QualifiedId__Group_1__0 )* ) )
            // InternalElkGraph.g:4436:1: ( ( rule__QualifiedId__Group_1__0 )* )
            {
            // InternalElkGraph.g:4436:1: ( ( rule__QualifiedId__Group_1__0 )* )
            // InternalElkGraph.g:4437:2: ( rule__QualifiedId__Group_1__0 )*
            {
             before(grammarAccess.getQualifiedIdAccess().getGroup_1()); 
            // InternalElkGraph.g:4438:2: ( rule__QualifiedId__Group_1__0 )*
            loop41:
            do {
                int alt41=2;
                int LA41_0 = input.LA(1);

                if ( (LA41_0==38) ) {
                    alt41=1;
                }


                switch (alt41) {
            	case 1 :
            	    // InternalElkGraph.g:4438:3: rule__QualifiedId__Group_1__0
            	    {
            	    pushFollow(FOLLOW_29);
            	    rule__QualifiedId__Group_1__0();

            	    state._fsp--;


            	    }
            	    break;

            	default :
            	    break loop41;
                }
            } while (true);

             after(grammarAccess.getQualifiedIdAccess().getGroup_1()); 

            }


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__QualifiedId__Group__1__Impl"


    // $ANTLR start "rule__QualifiedId__Group_1__0"
    // InternalElkGraph.g:4447:1: rule__QualifiedId__Group_1__0 : rule__QualifiedId__Group_1__0__Impl rule__QualifiedId__Group_1__1 ;
    public final void rule__QualifiedId__Group_1__0() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:4451:1: ( rule__QualifiedId__Group_1__0__Impl rule__QualifiedId__Group_1__1 )
            // InternalElkGraph.g:4452:2: rule__QualifiedId__Group_1__0__Impl rule__QualifiedId__Group_1__1
            {
            pushFollow(FOLLOW_7);
            rule__QualifiedId__Group_1__0__Impl();

            state._fsp--;

            pushFollow(FOLLOW_2);
            rule__QualifiedId__Group_1__1();

            state._fsp--;


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__QualifiedId__Group_1__0"


    // $ANTLR start "rule__QualifiedId__Group_1__0__Impl"
    // InternalElkGraph.g:4459:1: rule__QualifiedId__Group_1__0__Impl : ( '.' ) ;
    public final void rule__QualifiedId__Group_1__0__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:4463:1: ( ( '.' ) )
            // InternalElkGraph.g:4464:1: ( '.' )
            {
            // InternalElkGraph.g:4464:1: ( '.' )
            // InternalElkGraph.g:4465:2: '.'
            {
             before(grammarAccess.getQualifiedIdAccess().getFullStopKeyword_1_0()); 
            match(input,38,FOLLOW_2); 
             after(grammarAccess.getQualifiedIdAccess().getFullStopKeyword_1_0()); 

            }


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__QualifiedId__Group_1__0__Impl"


    // $ANTLR start "rule__QualifiedId__Group_1__1"
    // InternalElkGraph.g:4474:1: rule__QualifiedId__Group_1__1 : rule__QualifiedId__Group_1__1__Impl ;
    public final void rule__QualifiedId__Group_1__1() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:4478:1: ( rule__QualifiedId__Group_1__1__Impl )
            // InternalElkGraph.g:4479:2: rule__QualifiedId__Group_1__1__Impl
            {
            pushFollow(FOLLOW_2);
            rule__QualifiedId__Group_1__1__Impl();

            state._fsp--;


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__QualifiedId__Group_1__1"


    // $ANTLR start "rule__QualifiedId__Group_1__1__Impl"
    // InternalElkGraph.g:4485:1: rule__QualifiedId__Group_1__1__Impl : ( RULE_ID ) ;
    public final void rule__QualifiedId__Group_1__1__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:4489:1: ( ( RULE_ID ) )
            // InternalElkGraph.g:4490:1: ( RULE_ID )
            {
            // InternalElkGraph.g:4490:1: ( RULE_ID )
            // InternalElkGraph.g:4491:2: RULE_ID
            {
             before(grammarAccess.getQualifiedIdAccess().getIDTerminalRuleCall_1_1()); 
            match(input,RULE_ID,FOLLOW_2); 
             after(grammarAccess.getQualifiedIdAccess().getIDTerminalRuleCall_1_1()); 

            }


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__QualifiedId__Group_1__1__Impl"


    // $ANTLR start "rule__Property__Group__0"
    // InternalElkGraph.g:4501:1: rule__Property__Group__0 : rule__Property__Group__0__Impl rule__Property__Group__1 ;
    public final void rule__Property__Group__0() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:4505:1: ( rule__Property__Group__0__Impl rule__Property__Group__1 )
            // InternalElkGraph.g:4506:2: rule__Property__Group__0__Impl rule__Property__Group__1
            {
            pushFollow(FOLLOW_11);
            rule__Property__Group__0__Impl();

            state._fsp--;

            pushFollow(FOLLOW_2);
            rule__Property__Group__1();

            state._fsp--;


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__Property__Group__0"


    // $ANTLR start "rule__Property__Group__0__Impl"
    // InternalElkGraph.g:4513:1: rule__Property__Group__0__Impl : ( ( rule__Property__KeyAssignment_0 ) ) ;
    public final void rule__Property__Group__0__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:4517:1: ( ( ( rule__Property__KeyAssignment_0 ) ) )
            // InternalElkGraph.g:4518:1: ( ( rule__Property__KeyAssignment_0 ) )
            {
            // InternalElkGraph.g:4518:1: ( ( rule__Property__KeyAssignment_0 ) )
            // InternalElkGraph.g:4519:2: ( rule__Property__KeyAssignment_0 )
            {
             before(grammarAccess.getPropertyAccess().getKeyAssignment_0()); 
            // InternalElkGraph.g:4520:2: ( rule__Property__KeyAssignment_0 )
            // InternalElkGraph.g:4520:3: rule__Property__KeyAssignment_0
            {
            pushFollow(FOLLOW_2);
            rule__Property__KeyAssignment_0();

            state._fsp--;


            }

             after(grammarAccess.getPropertyAccess().getKeyAssignment_0()); 

            }


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__Property__Group__0__Impl"


    // $ANTLR start "rule__Property__Group__1"
    // InternalElkGraph.g:4528:1: rule__Property__Group__1 : rule__Property__Group__1__Impl rule__Property__Group__2 ;
    public final void rule__Property__Group__1() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:4532:1: ( rule__Property__Group__1__Impl rule__Property__Group__2 )
            // InternalElkGraph.g:4533:2: rule__Property__Group__1__Impl rule__Property__Group__2
            {
            pushFollow(FOLLOW_30);
            rule__Property__Group__1__Impl();

            state._fsp--;

            pushFollow(FOLLOW_2);
            rule__Property__Group__2();

            state._fsp--;


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__Property__Group__1"


    // $ANTLR start "rule__Property__Group__1__Impl"
    // InternalElkGraph.g:4540:1: rule__Property__Group__1__Impl : ( ':' ) ;
    public final void rule__Property__Group__1__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:4544:1: ( ( ':' ) )
            // InternalElkGraph.g:4545:1: ( ':' )
            {
            // InternalElkGraph.g:4545:1: ( ':' )
            // InternalElkGraph.g:4546:2: ':'
            {
             before(grammarAccess.getPropertyAccess().getColonKeyword_1()); 
            match(input,21,FOLLOW_2); 
             after(grammarAccess.getPropertyAccess().getColonKeyword_1()); 

            }


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__Property__Group__1__Impl"


    // $ANTLR start "rule__Property__Group__2"
    // InternalElkGraph.g:4555:1: rule__Property__Group__2 : rule__Property__Group__2__Impl ;
    public final void rule__Property__Group__2() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:4559:1: ( rule__Property__Group__2__Impl )
            // InternalElkGraph.g:4560:2: rule__Property__Group__2__Impl
            {
            pushFollow(FOLLOW_2);
            rule__Property__Group__2__Impl();

            state._fsp--;


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__Property__Group__2"


    // $ANTLR start "rule__Property__Group__2__Impl"
    // InternalElkGraph.g:4566:1: rule__Property__Group__2__Impl : ( ( rule__Property__Alternatives_2 ) ) ;
    public final void rule__Property__Group__2__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:4570:1: ( ( ( rule__Property__Alternatives_2 ) ) )
            // InternalElkGraph.g:4571:1: ( ( rule__Property__Alternatives_2 ) )
            {
            // InternalElkGraph.g:4571:1: ( ( rule__Property__Alternatives_2 ) )
            // InternalElkGraph.g:4572:2: ( rule__Property__Alternatives_2 )
            {
             before(grammarAccess.getPropertyAccess().getAlternatives_2()); 
            // InternalElkGraph.g:4573:2: ( rule__Property__Alternatives_2 )
            // InternalElkGraph.g:4573:3: rule__Property__Alternatives_2
            {
            pushFollow(FOLLOW_2);
            rule__Property__Alternatives_2();

            state._fsp--;


            }

             after(grammarAccess.getPropertyAccess().getAlternatives_2()); 

            }


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__Property__Group__2__Impl"


    // $ANTLR start "rule__PropertyKey__Group__0"
    // InternalElkGraph.g:4582:1: rule__PropertyKey__Group__0 : rule__PropertyKey__Group__0__Impl rule__PropertyKey__Group__1 ;
    public final void rule__PropertyKey__Group__0() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:4586:1: ( rule__PropertyKey__Group__0__Impl rule__PropertyKey__Group__1 )
            // InternalElkGraph.g:4587:2: rule__PropertyKey__Group__0__Impl rule__PropertyKey__Group__1
            {
            pushFollow(FOLLOW_28);
            rule__PropertyKey__Group__0__Impl();

            state._fsp--;

            pushFollow(FOLLOW_2);
            rule__PropertyKey__Group__1();

            state._fsp--;


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__PropertyKey__Group__0"


    // $ANTLR start "rule__PropertyKey__Group__0__Impl"
    // InternalElkGraph.g:4594:1: rule__PropertyKey__Group__0__Impl : ( RULE_ID ) ;
    public final void rule__PropertyKey__Group__0__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:4598:1: ( ( RULE_ID ) )
            // InternalElkGraph.g:4599:1: ( RULE_ID )
            {
            // InternalElkGraph.g:4599:1: ( RULE_ID )
            // InternalElkGraph.g:4600:2: RULE_ID
            {
             before(grammarAccess.getPropertyKeyAccess().getIDTerminalRuleCall_0()); 
            match(input,RULE_ID,FOLLOW_2); 
             after(grammarAccess.getPropertyKeyAccess().getIDTerminalRuleCall_0()); 

            }


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__PropertyKey__Group__0__Impl"


    // $ANTLR start "rule__PropertyKey__Group__1"
    // InternalElkGraph.g:4609:1: rule__PropertyKey__Group__1 : rule__PropertyKey__Group__1__Impl ;
    public final void rule__PropertyKey__Group__1() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:4613:1: ( rule__PropertyKey__Group__1__Impl )
            // InternalElkGraph.g:4614:2: rule__PropertyKey__Group__1__Impl
            {
            pushFollow(FOLLOW_2);
            rule__PropertyKey__Group__1__Impl();

            state._fsp--;


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__PropertyKey__Group__1"


    // $ANTLR start "rule__PropertyKey__Group__1__Impl"
    // InternalElkGraph.g:4620:1: rule__PropertyKey__Group__1__Impl : ( ( rule__PropertyKey__Group_1__0 )* ) ;
    public final void rule__PropertyKey__Group__1__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:4624:1: ( ( ( rule__PropertyKey__Group_1__0 )* ) )
            // InternalElkGraph.g:4625:1: ( ( rule__PropertyKey__Group_1__0 )* )
            {
            // InternalElkGraph.g:4625:1: ( ( rule__PropertyKey__Group_1__0 )* )
            // InternalElkGraph.g:4626:2: ( rule__PropertyKey__Group_1__0 )*
            {
             before(grammarAccess.getPropertyKeyAccess().getGroup_1()); 
            // InternalElkGraph.g:4627:2: ( rule__PropertyKey__Group_1__0 )*
            loop42:
            do {
                int alt42=2;
                int LA42_0 = input.LA(1);

                if ( (LA42_0==38) ) {
                    alt42=1;
                }


                switch (alt42) {
            	case 1 :
            	    // InternalElkGraph.g:4627:3: rule__PropertyKey__Group_1__0
            	    {
            	    pushFollow(FOLLOW_29);
            	    rule__PropertyKey__Group_1__0();

            	    state._fsp--;


            	    }
            	    break;

            	default :
            	    break loop42;
                }
            } while (true);

             after(grammarAccess.getPropertyKeyAccess().getGroup_1()); 

            }


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__PropertyKey__Group__1__Impl"


    // $ANTLR start "rule__PropertyKey__Group_1__0"
    // InternalElkGraph.g:4636:1: rule__PropertyKey__Group_1__0 : rule__PropertyKey__Group_1__0__Impl rule__PropertyKey__Group_1__1 ;
    public final void rule__PropertyKey__Group_1__0() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:4640:1: ( rule__PropertyKey__Group_1__0__Impl rule__PropertyKey__Group_1__1 )
            // InternalElkGraph.g:4641:2: rule__PropertyKey__Group_1__0__Impl rule__PropertyKey__Group_1__1
            {
            pushFollow(FOLLOW_7);
            rule__PropertyKey__Group_1__0__Impl();

            state._fsp--;

            pushFollow(FOLLOW_2);
            rule__PropertyKey__Group_1__1();

            state._fsp--;


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__PropertyKey__Group_1__0"


    // $ANTLR start "rule__PropertyKey__Group_1__0__Impl"
    // InternalElkGraph.g:4648:1: rule__PropertyKey__Group_1__0__Impl : ( '.' ) ;
    public final void rule__PropertyKey__Group_1__0__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:4652:1: ( ( '.' ) )
            // InternalElkGraph.g:4653:1: ( '.' )
            {
            // InternalElkGraph.g:4653:1: ( '.' )
            // InternalElkGraph.g:4654:2: '.'
            {
             before(grammarAccess.getPropertyKeyAccess().getFullStopKeyword_1_0()); 
            match(input,38,FOLLOW_2); 
             after(grammarAccess.getPropertyKeyAccess().getFullStopKeyword_1_0()); 

            }


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__PropertyKey__Group_1__0__Impl"


    // $ANTLR start "rule__PropertyKey__Group_1__1"
    // InternalElkGraph.g:4663:1: rule__PropertyKey__Group_1__1 : rule__PropertyKey__Group_1__1__Impl ;
    public final void rule__PropertyKey__Group_1__1() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:4667:1: ( rule__PropertyKey__Group_1__1__Impl )
            // InternalElkGraph.g:4668:2: rule__PropertyKey__Group_1__1__Impl
            {
            pushFollow(FOLLOW_2);
            rule__PropertyKey__Group_1__1__Impl();

            state._fsp--;


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__PropertyKey__Group_1__1"


    // $ANTLR start "rule__PropertyKey__Group_1__1__Impl"
    // InternalElkGraph.g:4674:1: rule__PropertyKey__Group_1__1__Impl : ( RULE_ID ) ;
    public final void rule__PropertyKey__Group_1__1__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:4678:1: ( ( RULE_ID ) )
            // InternalElkGraph.g:4679:1: ( RULE_ID )
            {
            // InternalElkGraph.g:4679:1: ( RULE_ID )
            // InternalElkGraph.g:4680:2: RULE_ID
            {
             before(grammarAccess.getPropertyKeyAccess().getIDTerminalRuleCall_1_1()); 
            match(input,RULE_ID,FOLLOW_2); 
             after(grammarAccess.getPropertyKeyAccess().getIDTerminalRuleCall_1_1()); 

            }


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__PropertyKey__Group_1__1__Impl"


    // $ANTLR start "rule__ShapeLayout__UnorderedGroup_2"
    // InternalElkGraph.g:4690:1: rule__ShapeLayout__UnorderedGroup_2 : ( rule__ShapeLayout__UnorderedGroup_2__0 )? ;
    public final void rule__ShapeLayout__UnorderedGroup_2() throws RecognitionException {

        		int stackSize = keepStackSize();
        		getUnorderedGroupHelper().enter(grammarAccess.getShapeLayoutAccess().getUnorderedGroup_2());
        	
        try {
            // InternalElkGraph.g:4695:1: ( ( rule__ShapeLayout__UnorderedGroup_2__0 )? )
            // InternalElkGraph.g:4696:2: ( rule__ShapeLayout__UnorderedGroup_2__0 )?
            {
            // InternalElkGraph.g:4696:2: ( rule__ShapeLayout__UnorderedGroup_2__0 )?
            int alt43=2;
            int LA43_0 = input.LA(1);

            if ( LA43_0 == 26 && getUnorderedGroupHelper().canSelect(grammarAccess.getShapeLayoutAccess().getUnorderedGroup_2(), 0) ) {
                alt43=1;
            }
            else if ( LA43_0 == 28 && getUnorderedGroupHelper().canSelect(grammarAccess.getShapeLayoutAccess().getUnorderedGroup_2(), 1) ) {
                alt43=1;
            }
            switch (alt43) {
                case 1 :
                    // InternalElkGraph.g:4696:2: rule__ShapeLayout__UnorderedGroup_2__0
                    {
                    pushFollow(FOLLOW_2);
                    rule__ShapeLayout__UnorderedGroup_2__0();

                    state._fsp--;


                    }
                    break;

            }


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	getUnorderedGroupHelper().leave(grammarAccess.getShapeLayoutAccess().getUnorderedGroup_2());
            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ShapeLayout__UnorderedGroup_2"


    // $ANTLR start "rule__ShapeLayout__UnorderedGroup_2__Impl"
    // InternalElkGraph.g:4704:1: rule__ShapeLayout__UnorderedGroup_2__Impl : ( ({...}? => ( ( ( rule__ShapeLayout__Group_2_0__0 ) ) ) ) | ({...}? => ( ( ( rule__ShapeLayout__Group_2_1__0 ) ) ) ) ) ;
    public final void rule__ShapeLayout__UnorderedGroup_2__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        		boolean selected = false;
        	
        try {
            // InternalElkGraph.g:4709:1: ( ( ({...}? => ( ( ( rule__ShapeLayout__Group_2_0__0 ) ) ) ) | ({...}? => ( ( ( rule__ShapeLayout__Group_2_1__0 ) ) ) ) ) )
            // InternalElkGraph.g:4710:3: ( ({...}? => ( ( ( rule__ShapeLayout__Group_2_0__0 ) ) ) ) | ({...}? => ( ( ( rule__ShapeLayout__Group_2_1__0 ) ) ) ) )
            {
            // InternalElkGraph.g:4710:3: ( ({...}? => ( ( ( rule__ShapeLayout__Group_2_0__0 ) ) ) ) | ({...}? => ( ( ( rule__ShapeLayout__Group_2_1__0 ) ) ) ) )
            int alt44=2;
            int LA44_0 = input.LA(1);

            if ( LA44_0 == 26 && getUnorderedGroupHelper().canSelect(grammarAccess.getShapeLayoutAccess().getUnorderedGroup_2(), 0) ) {
                alt44=1;
            }
            else if ( LA44_0 == 28 && getUnorderedGroupHelper().canSelect(grammarAccess.getShapeLayoutAccess().getUnorderedGroup_2(), 1) ) {
                alt44=2;
            }
            else {
                NoViableAltException nvae =
                    new NoViableAltException("", 44, 0, input);

                throw nvae;
            }
            switch (alt44) {
                case 1 :
                    // InternalElkGraph.g:4711:3: ({...}? => ( ( ( rule__ShapeLayout__Group_2_0__0 ) ) ) )
                    {
                    // InternalElkGraph.g:4711:3: ({...}? => ( ( ( rule__ShapeLayout__Group_2_0__0 ) ) ) )
                    // InternalElkGraph.g:4712:4: {...}? => ( ( ( rule__ShapeLayout__Group_2_0__0 ) ) )
                    {
                    if ( ! getUnorderedGroupHelper().canSelect(grammarAccess.getShapeLayoutAccess().getUnorderedGroup_2(), 0) ) {
                        throw new FailedPredicateException(input, "rule__ShapeLayout__UnorderedGroup_2__Impl", "getUnorderedGroupHelper().canSelect(grammarAccess.getShapeLayoutAccess().getUnorderedGroup_2(), 0)");
                    }
                    // InternalElkGraph.g:4712:107: ( ( ( rule__ShapeLayout__Group_2_0__0 ) ) )
                    // InternalElkGraph.g:4713:5: ( ( rule__ShapeLayout__Group_2_0__0 ) )
                    {

                    					getUnorderedGroupHelper().select(grammarAccess.getShapeLayoutAccess().getUnorderedGroup_2(), 0);
                    				

                    					selected = true;
                    				
                    // InternalElkGraph.g:4719:5: ( ( rule__ShapeLayout__Group_2_0__0 ) )
                    // InternalElkGraph.g:4720:6: ( rule__ShapeLayout__Group_2_0__0 )
                    {
                     before(grammarAccess.getShapeLayoutAccess().getGroup_2_0()); 
                    // InternalElkGraph.g:4721:6: ( rule__ShapeLayout__Group_2_0__0 )
                    // InternalElkGraph.g:4721:7: rule__ShapeLayout__Group_2_0__0
                    {
                    pushFollow(FOLLOW_2);
                    rule__ShapeLayout__Group_2_0__0();

                    state._fsp--;


                    }

                     after(grammarAccess.getShapeLayoutAccess().getGroup_2_0()); 

                    }


                    }


                    }


                    }
                    break;
                case 2 :
                    // InternalElkGraph.g:4726:3: ({...}? => ( ( ( rule__ShapeLayout__Group_2_1__0 ) ) ) )
                    {
                    // InternalElkGraph.g:4726:3: ({...}? => ( ( ( rule__ShapeLayout__Group_2_1__0 ) ) ) )
                    // InternalElkGraph.g:4727:4: {...}? => ( ( ( rule__ShapeLayout__Group_2_1__0 ) ) )
                    {
                    if ( ! getUnorderedGroupHelper().canSelect(grammarAccess.getShapeLayoutAccess().getUnorderedGroup_2(), 1) ) {
                        throw new FailedPredicateException(input, "rule__ShapeLayout__UnorderedGroup_2__Impl", "getUnorderedGroupHelper().canSelect(grammarAccess.getShapeLayoutAccess().getUnorderedGroup_2(), 1)");
                    }
                    // InternalElkGraph.g:4727:107: ( ( ( rule__ShapeLayout__Group_2_1__0 ) ) )
                    // InternalElkGraph.g:4728:5: ( ( rule__ShapeLayout__Group_2_1__0 ) )
                    {

                    					getUnorderedGroupHelper().select(grammarAccess.getShapeLayoutAccess().getUnorderedGroup_2(), 1);
                    				

                    					selected = true;
                    				
                    // InternalElkGraph.g:4734:5: ( ( rule__ShapeLayout__Group_2_1__0 ) )
                    // InternalElkGraph.g:4735:6: ( rule__ShapeLayout__Group_2_1__0 )
                    {
                     before(grammarAccess.getShapeLayoutAccess().getGroup_2_1()); 
                    // InternalElkGraph.g:4736:6: ( rule__ShapeLayout__Group_2_1__0 )
                    // InternalElkGraph.g:4736:7: rule__ShapeLayout__Group_2_1__0
                    {
                    pushFollow(FOLLOW_2);
                    rule__ShapeLayout__Group_2_1__0();

                    state._fsp--;


                    }

                     after(grammarAccess.getShapeLayoutAccess().getGroup_2_1()); 

                    }


                    }


                    }


                    }
                    break;

            }


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	if (selected)
            		getUnorderedGroupHelper().returnFromSelection(grammarAccess.getShapeLayoutAccess().getUnorderedGroup_2());
            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ShapeLayout__UnorderedGroup_2__Impl"


    // $ANTLR start "rule__ShapeLayout__UnorderedGroup_2__0"
    // InternalElkGraph.g:4749:1: rule__ShapeLayout__UnorderedGroup_2__0 : rule__ShapeLayout__UnorderedGroup_2__Impl ( rule__ShapeLayout__UnorderedGroup_2__1 )? ;
    public final void rule__ShapeLayout__UnorderedGroup_2__0() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:4753:1: ( rule__ShapeLayout__UnorderedGroup_2__Impl ( rule__ShapeLayout__UnorderedGroup_2__1 )? )
            // InternalElkGraph.g:4754:2: rule__ShapeLayout__UnorderedGroup_2__Impl ( rule__ShapeLayout__UnorderedGroup_2__1 )?
            {
            pushFollow(FOLLOW_31);
            rule__ShapeLayout__UnorderedGroup_2__Impl();

            state._fsp--;

            // InternalElkGraph.g:4755:2: ( rule__ShapeLayout__UnorderedGroup_2__1 )?
            int alt45=2;
            int LA45_0 = input.LA(1);

            if ( LA45_0 == 26 && getUnorderedGroupHelper().canSelect(grammarAccess.getShapeLayoutAccess().getUnorderedGroup_2(), 0) ) {
                alt45=1;
            }
            else if ( LA45_0 == 28 && getUnorderedGroupHelper().canSelect(grammarAccess.getShapeLayoutAccess().getUnorderedGroup_2(), 1) ) {
                alt45=1;
            }
            switch (alt45) {
                case 1 :
                    // InternalElkGraph.g:4755:2: rule__ShapeLayout__UnorderedGroup_2__1
                    {
                    pushFollow(FOLLOW_2);
                    rule__ShapeLayout__UnorderedGroup_2__1();

                    state._fsp--;


                    }
                    break;

            }


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ShapeLayout__UnorderedGroup_2__0"


    // $ANTLR start "rule__ShapeLayout__UnorderedGroup_2__1"
    // InternalElkGraph.g:4761:1: rule__ShapeLayout__UnorderedGroup_2__1 : rule__ShapeLayout__UnorderedGroup_2__Impl ;
    public final void rule__ShapeLayout__UnorderedGroup_2__1() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:4765:1: ( rule__ShapeLayout__UnorderedGroup_2__Impl )
            // InternalElkGraph.g:4766:2: rule__ShapeLayout__UnorderedGroup_2__Impl
            {
            pushFollow(FOLLOW_2);
            rule__ShapeLayout__UnorderedGroup_2__Impl();

            state._fsp--;


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ShapeLayout__UnorderedGroup_2__1"


    // $ANTLR start "rule__ElkSingleEdgeSection__UnorderedGroup_1_0"
    // InternalElkGraph.g:4773:1: rule__ElkSingleEdgeSection__UnorderedGroup_1_0 : ( rule__ElkSingleEdgeSection__UnorderedGroup_1_0__0 )? ;
    public final void rule__ElkSingleEdgeSection__UnorderedGroup_1_0() throws RecognitionException {

        		int stackSize = keepStackSize();
        		getUnorderedGroupHelper().enter(grammarAccess.getElkSingleEdgeSectionAccess().getUnorderedGroup_1_0());
        	
        try {
            // InternalElkGraph.g:4778:1: ( ( rule__ElkSingleEdgeSection__UnorderedGroup_1_0__0 )? )
            // InternalElkGraph.g:4779:2: ( rule__ElkSingleEdgeSection__UnorderedGroup_1_0__0 )?
            {
            // InternalElkGraph.g:4779:2: ( rule__ElkSingleEdgeSection__UnorderedGroup_1_0__0 )?
            int alt46=2;
            int LA46_0 = input.LA(1);

            if ( LA46_0 == 31 && getUnorderedGroupHelper().canSelect(grammarAccess.getElkSingleEdgeSectionAccess().getUnorderedGroup_1_0(), 0) ) {
                alt46=1;
            }
            else if ( LA46_0 == 32 && getUnorderedGroupHelper().canSelect(grammarAccess.getElkSingleEdgeSectionAccess().getUnorderedGroup_1_0(), 1) ) {
                alt46=1;
            }
            else if ( LA46_0 == 33 && getUnorderedGroupHelper().canSelect(grammarAccess.getElkSingleEdgeSectionAccess().getUnorderedGroup_1_0(), 2) ) {
                alt46=1;
            }
            else if ( LA46_0 == 34 && getUnorderedGroupHelper().canSelect(grammarAccess.getElkSingleEdgeSectionAccess().getUnorderedGroup_1_0(), 3) ) {
                alt46=1;
            }
            switch (alt46) {
                case 1 :
                    // InternalElkGraph.g:4779:2: rule__ElkSingleEdgeSection__UnorderedGroup_1_0__0
                    {
                    pushFollow(FOLLOW_2);
                    rule__ElkSingleEdgeSection__UnorderedGroup_1_0__0();

                    state._fsp--;


                    }
                    break;

            }


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	getUnorderedGroupHelper().leave(grammarAccess.getElkSingleEdgeSectionAccess().getUnorderedGroup_1_0());
            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ElkSingleEdgeSection__UnorderedGroup_1_0"


    // $ANTLR start "rule__ElkSingleEdgeSection__UnorderedGroup_1_0__Impl"
    // InternalElkGraph.g:4787:1: rule__ElkSingleEdgeSection__UnorderedGroup_1_0__Impl : ( ({...}? => ( ( ( rule__ElkSingleEdgeSection__Group_1_0_0__0 ) ) ) ) | ({...}? => ( ( ( rule__ElkSingleEdgeSection__Group_1_0_1__0 ) ) ) ) | ({...}? => ( ( ( rule__ElkSingleEdgeSection__Group_1_0_2__0 ) ) ) ) | ({...}? => ( ( ( rule__ElkSingleEdgeSection__Group_1_0_3__0 ) ) ) ) ) ;
    public final void rule__ElkSingleEdgeSection__UnorderedGroup_1_0__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        		boolean selected = false;
        	
        try {
            // InternalElkGraph.g:4792:1: ( ( ({...}? => ( ( ( rule__ElkSingleEdgeSection__Group_1_0_0__0 ) ) ) ) | ({...}? => ( ( ( rule__ElkSingleEdgeSection__Group_1_0_1__0 ) ) ) ) | ({...}? => ( ( ( rule__ElkSingleEdgeSection__Group_1_0_2__0 ) ) ) ) | ({...}? => ( ( ( rule__ElkSingleEdgeSection__Group_1_0_3__0 ) ) ) ) ) )
            // InternalElkGraph.g:4793:3: ( ({...}? => ( ( ( rule__ElkSingleEdgeSection__Group_1_0_0__0 ) ) ) ) | ({...}? => ( ( ( rule__ElkSingleEdgeSection__Group_1_0_1__0 ) ) ) ) | ({...}? => ( ( ( rule__ElkSingleEdgeSection__Group_1_0_2__0 ) ) ) ) | ({...}? => ( ( ( rule__ElkSingleEdgeSection__Group_1_0_3__0 ) ) ) ) )
            {
            // InternalElkGraph.g:4793:3: ( ({...}? => ( ( ( rule__ElkSingleEdgeSection__Group_1_0_0__0 ) ) ) ) | ({...}? => ( ( ( rule__ElkSingleEdgeSection__Group_1_0_1__0 ) ) ) ) | ({...}? => ( ( ( rule__ElkSingleEdgeSection__Group_1_0_2__0 ) ) ) ) | ({...}? => ( ( ( rule__ElkSingleEdgeSection__Group_1_0_3__0 ) ) ) ) )
            int alt47=4;
            int LA47_0 = input.LA(1);

            if ( LA47_0 == 31 && getUnorderedGroupHelper().canSelect(grammarAccess.getElkSingleEdgeSectionAccess().getUnorderedGroup_1_0(), 0) ) {
                alt47=1;
            }
            else if ( LA47_0 == 32 && getUnorderedGroupHelper().canSelect(grammarAccess.getElkSingleEdgeSectionAccess().getUnorderedGroup_1_0(), 1) ) {
                alt47=2;
            }
            else if ( LA47_0 == 33 && getUnorderedGroupHelper().canSelect(grammarAccess.getElkSingleEdgeSectionAccess().getUnorderedGroup_1_0(), 2) ) {
                alt47=3;
            }
            else if ( LA47_0 == 34 && getUnorderedGroupHelper().canSelect(grammarAccess.getElkSingleEdgeSectionAccess().getUnorderedGroup_1_0(), 3) ) {
                alt47=4;
            }
            else {
                NoViableAltException nvae =
                    new NoViableAltException("", 47, 0, input);

                throw nvae;
            }
            switch (alt47) {
                case 1 :
                    // InternalElkGraph.g:4794:3: ({...}? => ( ( ( rule__ElkSingleEdgeSection__Group_1_0_0__0 ) ) ) )
                    {
                    // InternalElkGraph.g:4794:3: ({...}? => ( ( ( rule__ElkSingleEdgeSection__Group_1_0_0__0 ) ) ) )
                    // InternalElkGraph.g:4795:4: {...}? => ( ( ( rule__ElkSingleEdgeSection__Group_1_0_0__0 ) ) )
                    {
                    if ( ! getUnorderedGroupHelper().canSelect(grammarAccess.getElkSingleEdgeSectionAccess().getUnorderedGroup_1_0(), 0) ) {
                        throw new FailedPredicateException(input, "rule__ElkSingleEdgeSection__UnorderedGroup_1_0__Impl", "getUnorderedGroupHelper().canSelect(grammarAccess.getElkSingleEdgeSectionAccess().getUnorderedGroup_1_0(), 0)");
                    }
                    // InternalElkGraph.g:4795:118: ( ( ( rule__ElkSingleEdgeSection__Group_1_0_0__0 ) ) )
                    // InternalElkGraph.g:4796:5: ( ( rule__ElkSingleEdgeSection__Group_1_0_0__0 ) )
                    {

                    					getUnorderedGroupHelper().select(grammarAccess.getElkSingleEdgeSectionAccess().getUnorderedGroup_1_0(), 0);
                    				

                    					selected = true;
                    				
                    // InternalElkGraph.g:4802:5: ( ( rule__ElkSingleEdgeSection__Group_1_0_0__0 ) )
                    // InternalElkGraph.g:4803:6: ( rule__ElkSingleEdgeSection__Group_1_0_0__0 )
                    {
                     before(grammarAccess.getElkSingleEdgeSectionAccess().getGroup_1_0_0()); 
                    // InternalElkGraph.g:4804:6: ( rule__ElkSingleEdgeSection__Group_1_0_0__0 )
                    // InternalElkGraph.g:4804:7: rule__ElkSingleEdgeSection__Group_1_0_0__0
                    {
                    pushFollow(FOLLOW_2);
                    rule__ElkSingleEdgeSection__Group_1_0_0__0();

                    state._fsp--;


                    }

                     after(grammarAccess.getElkSingleEdgeSectionAccess().getGroup_1_0_0()); 

                    }


                    }


                    }


                    }
                    break;
                case 2 :
                    // InternalElkGraph.g:4809:3: ({...}? => ( ( ( rule__ElkSingleEdgeSection__Group_1_0_1__0 ) ) ) )
                    {
                    // InternalElkGraph.g:4809:3: ({...}? => ( ( ( rule__ElkSingleEdgeSection__Group_1_0_1__0 ) ) ) )
                    // InternalElkGraph.g:4810:4: {...}? => ( ( ( rule__ElkSingleEdgeSection__Group_1_0_1__0 ) ) )
                    {
                    if ( ! getUnorderedGroupHelper().canSelect(grammarAccess.getElkSingleEdgeSectionAccess().getUnorderedGroup_1_0(), 1) ) {
                        throw new FailedPredicateException(input, "rule__ElkSingleEdgeSection__UnorderedGroup_1_0__Impl", "getUnorderedGroupHelper().canSelect(grammarAccess.getElkSingleEdgeSectionAccess().getUnorderedGroup_1_0(), 1)");
                    }
                    // InternalElkGraph.g:4810:118: ( ( ( rule__ElkSingleEdgeSection__Group_1_0_1__0 ) ) )
                    // InternalElkGraph.g:4811:5: ( ( rule__ElkSingleEdgeSection__Group_1_0_1__0 ) )
                    {

                    					getUnorderedGroupHelper().select(grammarAccess.getElkSingleEdgeSectionAccess().getUnorderedGroup_1_0(), 1);
                    				

                    					selected = true;
                    				
                    // InternalElkGraph.g:4817:5: ( ( rule__ElkSingleEdgeSection__Group_1_0_1__0 ) )
                    // InternalElkGraph.g:4818:6: ( rule__ElkSingleEdgeSection__Group_1_0_1__0 )
                    {
                     before(grammarAccess.getElkSingleEdgeSectionAccess().getGroup_1_0_1()); 
                    // InternalElkGraph.g:4819:6: ( rule__ElkSingleEdgeSection__Group_1_0_1__0 )
                    // InternalElkGraph.g:4819:7: rule__ElkSingleEdgeSection__Group_1_0_1__0
                    {
                    pushFollow(FOLLOW_2);
                    rule__ElkSingleEdgeSection__Group_1_0_1__0();

                    state._fsp--;


                    }

                     after(grammarAccess.getElkSingleEdgeSectionAccess().getGroup_1_0_1()); 

                    }


                    }


                    }


                    }
                    break;
                case 3 :
                    // InternalElkGraph.g:4824:3: ({...}? => ( ( ( rule__ElkSingleEdgeSection__Group_1_0_2__0 ) ) ) )
                    {
                    // InternalElkGraph.g:4824:3: ({...}? => ( ( ( rule__ElkSingleEdgeSection__Group_1_0_2__0 ) ) ) )
                    // InternalElkGraph.g:4825:4: {...}? => ( ( ( rule__ElkSingleEdgeSection__Group_1_0_2__0 ) ) )
                    {
                    if ( ! getUnorderedGroupHelper().canSelect(grammarAccess.getElkSingleEdgeSectionAccess().getUnorderedGroup_1_0(), 2) ) {
                        throw new FailedPredicateException(input, "rule__ElkSingleEdgeSection__UnorderedGroup_1_0__Impl", "getUnorderedGroupHelper().canSelect(grammarAccess.getElkSingleEdgeSectionAccess().getUnorderedGroup_1_0(), 2)");
                    }
                    // InternalElkGraph.g:4825:118: ( ( ( rule__ElkSingleEdgeSection__Group_1_0_2__0 ) ) )
                    // InternalElkGraph.g:4826:5: ( ( rule__ElkSingleEdgeSection__Group_1_0_2__0 ) )
                    {

                    					getUnorderedGroupHelper().select(grammarAccess.getElkSingleEdgeSectionAccess().getUnorderedGroup_1_0(), 2);
                    				

                    					selected = true;
                    				
                    // InternalElkGraph.g:4832:5: ( ( rule__ElkSingleEdgeSection__Group_1_0_2__0 ) )
                    // InternalElkGraph.g:4833:6: ( rule__ElkSingleEdgeSection__Group_1_0_2__0 )
                    {
                     before(grammarAccess.getElkSingleEdgeSectionAccess().getGroup_1_0_2()); 
                    // InternalElkGraph.g:4834:6: ( rule__ElkSingleEdgeSection__Group_1_0_2__0 )
                    // InternalElkGraph.g:4834:7: rule__ElkSingleEdgeSection__Group_1_0_2__0
                    {
                    pushFollow(FOLLOW_2);
                    rule__ElkSingleEdgeSection__Group_1_0_2__0();

                    state._fsp--;


                    }

                     after(grammarAccess.getElkSingleEdgeSectionAccess().getGroup_1_0_2()); 

                    }


                    }


                    }


                    }
                    break;
                case 4 :
                    // InternalElkGraph.g:4839:3: ({...}? => ( ( ( rule__ElkSingleEdgeSection__Group_1_0_3__0 ) ) ) )
                    {
                    // InternalElkGraph.g:4839:3: ({...}? => ( ( ( rule__ElkSingleEdgeSection__Group_1_0_3__0 ) ) ) )
                    // InternalElkGraph.g:4840:4: {...}? => ( ( ( rule__ElkSingleEdgeSection__Group_1_0_3__0 ) ) )
                    {
                    if ( ! getUnorderedGroupHelper().canSelect(grammarAccess.getElkSingleEdgeSectionAccess().getUnorderedGroup_1_0(), 3) ) {
                        throw new FailedPredicateException(input, "rule__ElkSingleEdgeSection__UnorderedGroup_1_0__Impl", "getUnorderedGroupHelper().canSelect(grammarAccess.getElkSingleEdgeSectionAccess().getUnorderedGroup_1_0(), 3)");
                    }
                    // InternalElkGraph.g:4840:118: ( ( ( rule__ElkSingleEdgeSection__Group_1_0_3__0 ) ) )
                    // InternalElkGraph.g:4841:5: ( ( rule__ElkSingleEdgeSection__Group_1_0_3__0 ) )
                    {

                    					getUnorderedGroupHelper().select(grammarAccess.getElkSingleEdgeSectionAccess().getUnorderedGroup_1_0(), 3);
                    				

                    					selected = true;
                    				
                    // InternalElkGraph.g:4847:5: ( ( rule__ElkSingleEdgeSection__Group_1_0_3__0 ) )
                    // InternalElkGraph.g:4848:6: ( rule__ElkSingleEdgeSection__Group_1_0_3__0 )
                    {
                     before(grammarAccess.getElkSingleEdgeSectionAccess().getGroup_1_0_3()); 
                    // InternalElkGraph.g:4849:6: ( rule__ElkSingleEdgeSection__Group_1_0_3__0 )
                    // InternalElkGraph.g:4849:7: rule__ElkSingleEdgeSection__Group_1_0_3__0
                    {
                    pushFollow(FOLLOW_2);
                    rule__ElkSingleEdgeSection__Group_1_0_3__0();

                    state._fsp--;


                    }

                     after(grammarAccess.getElkSingleEdgeSectionAccess().getGroup_1_0_3()); 

                    }


                    }


                    }


                    }
                    break;

            }


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	if (selected)
            		getUnorderedGroupHelper().returnFromSelection(grammarAccess.getElkSingleEdgeSectionAccess().getUnorderedGroup_1_0());
            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ElkSingleEdgeSection__UnorderedGroup_1_0__Impl"


    // $ANTLR start "rule__ElkSingleEdgeSection__UnorderedGroup_1_0__0"
    // InternalElkGraph.g:4862:1: rule__ElkSingleEdgeSection__UnorderedGroup_1_0__0 : rule__ElkSingleEdgeSection__UnorderedGroup_1_0__Impl ( rule__ElkSingleEdgeSection__UnorderedGroup_1_0__1 )? ;
    public final void rule__ElkSingleEdgeSection__UnorderedGroup_1_0__0() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:4866:1: ( rule__ElkSingleEdgeSection__UnorderedGroup_1_0__Impl ( rule__ElkSingleEdgeSection__UnorderedGroup_1_0__1 )? )
            // InternalElkGraph.g:4867:2: rule__ElkSingleEdgeSection__UnorderedGroup_1_0__Impl ( rule__ElkSingleEdgeSection__UnorderedGroup_1_0__1 )?
            {
            pushFollow(FOLLOW_32);
            rule__ElkSingleEdgeSection__UnorderedGroup_1_0__Impl();

            state._fsp--;

            // InternalElkGraph.g:4868:2: ( rule__ElkSingleEdgeSection__UnorderedGroup_1_0__1 )?
            int alt48=2;
            int LA48_0 = input.LA(1);

            if ( LA48_0 == 31 && getUnorderedGroupHelper().canSelect(grammarAccess.getElkSingleEdgeSectionAccess().getUnorderedGroup_1_0(), 0) ) {
                alt48=1;
            }
            else if ( LA48_0 == 32 && getUnorderedGroupHelper().canSelect(grammarAccess.getElkSingleEdgeSectionAccess().getUnorderedGroup_1_0(), 1) ) {
                alt48=1;
            }
            else if ( LA48_0 == 33 && getUnorderedGroupHelper().canSelect(grammarAccess.getElkSingleEdgeSectionAccess().getUnorderedGroup_1_0(), 2) ) {
                alt48=1;
            }
            else if ( LA48_0 == 34 && getUnorderedGroupHelper().canSelect(grammarAccess.getElkSingleEdgeSectionAccess().getUnorderedGroup_1_0(), 3) ) {
                alt48=1;
            }
            switch (alt48) {
                case 1 :
                    // InternalElkGraph.g:4868:2: rule__ElkSingleEdgeSection__UnorderedGroup_1_0__1
                    {
                    pushFollow(FOLLOW_2);
                    rule__ElkSingleEdgeSection__UnorderedGroup_1_0__1();

                    state._fsp--;


                    }
                    break;

            }


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ElkSingleEdgeSection__UnorderedGroup_1_0__0"


    // $ANTLR start "rule__ElkSingleEdgeSection__UnorderedGroup_1_0__1"
    // InternalElkGraph.g:4874:1: rule__ElkSingleEdgeSection__UnorderedGroup_1_0__1 : rule__ElkSingleEdgeSection__UnorderedGroup_1_0__Impl ( rule__ElkSingleEdgeSection__UnorderedGroup_1_0__2 )? ;
    public final void rule__ElkSingleEdgeSection__UnorderedGroup_1_0__1() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:4878:1: ( rule__ElkSingleEdgeSection__UnorderedGroup_1_0__Impl ( rule__ElkSingleEdgeSection__UnorderedGroup_1_0__2 )? )
            // InternalElkGraph.g:4879:2: rule__ElkSingleEdgeSection__UnorderedGroup_1_0__Impl ( rule__ElkSingleEdgeSection__UnorderedGroup_1_0__2 )?
            {
            pushFollow(FOLLOW_32);
            rule__ElkSingleEdgeSection__UnorderedGroup_1_0__Impl();

            state._fsp--;

            // InternalElkGraph.g:4880:2: ( rule__ElkSingleEdgeSection__UnorderedGroup_1_0__2 )?
            int alt49=2;
            int LA49_0 = input.LA(1);

            if ( LA49_0 == 31 && getUnorderedGroupHelper().canSelect(grammarAccess.getElkSingleEdgeSectionAccess().getUnorderedGroup_1_0(), 0) ) {
                alt49=1;
            }
            else if ( LA49_0 == 32 && getUnorderedGroupHelper().canSelect(grammarAccess.getElkSingleEdgeSectionAccess().getUnorderedGroup_1_0(), 1) ) {
                alt49=1;
            }
            else if ( LA49_0 == 33 && getUnorderedGroupHelper().canSelect(grammarAccess.getElkSingleEdgeSectionAccess().getUnorderedGroup_1_0(), 2) ) {
                alt49=1;
            }
            else if ( LA49_0 == 34 && getUnorderedGroupHelper().canSelect(grammarAccess.getElkSingleEdgeSectionAccess().getUnorderedGroup_1_0(), 3) ) {
                alt49=1;
            }
            switch (alt49) {
                case 1 :
                    // InternalElkGraph.g:4880:2: rule__ElkSingleEdgeSection__UnorderedGroup_1_0__2
                    {
                    pushFollow(FOLLOW_2);
                    rule__ElkSingleEdgeSection__UnorderedGroup_1_0__2();

                    state._fsp--;


                    }
                    break;

            }


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ElkSingleEdgeSection__UnorderedGroup_1_0__1"


    // $ANTLR start "rule__ElkSingleEdgeSection__UnorderedGroup_1_0__2"
    // InternalElkGraph.g:4886:1: rule__ElkSingleEdgeSection__UnorderedGroup_1_0__2 : rule__ElkSingleEdgeSection__UnorderedGroup_1_0__Impl ( rule__ElkSingleEdgeSection__UnorderedGroup_1_0__3 )? ;
    public final void rule__ElkSingleEdgeSection__UnorderedGroup_1_0__2() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:4890:1: ( rule__ElkSingleEdgeSection__UnorderedGroup_1_0__Impl ( rule__ElkSingleEdgeSection__UnorderedGroup_1_0__3 )? )
            // InternalElkGraph.g:4891:2: rule__ElkSingleEdgeSection__UnorderedGroup_1_0__Impl ( rule__ElkSingleEdgeSection__UnorderedGroup_1_0__3 )?
            {
            pushFollow(FOLLOW_32);
            rule__ElkSingleEdgeSection__UnorderedGroup_1_0__Impl();

            state._fsp--;

            // InternalElkGraph.g:4892:2: ( rule__ElkSingleEdgeSection__UnorderedGroup_1_0__3 )?
            int alt50=2;
            int LA50_0 = input.LA(1);

            if ( LA50_0 == 31 && getUnorderedGroupHelper().canSelect(grammarAccess.getElkSingleEdgeSectionAccess().getUnorderedGroup_1_0(), 0) ) {
                alt50=1;
            }
            else if ( LA50_0 == 32 && getUnorderedGroupHelper().canSelect(grammarAccess.getElkSingleEdgeSectionAccess().getUnorderedGroup_1_0(), 1) ) {
                alt50=1;
            }
            else if ( LA50_0 == 33 && getUnorderedGroupHelper().canSelect(grammarAccess.getElkSingleEdgeSectionAccess().getUnorderedGroup_1_0(), 2) ) {
                alt50=1;
            }
            else if ( LA50_0 == 34 && getUnorderedGroupHelper().canSelect(grammarAccess.getElkSingleEdgeSectionAccess().getUnorderedGroup_1_0(), 3) ) {
                alt50=1;
            }
            switch (alt50) {
                case 1 :
                    // InternalElkGraph.g:4892:2: rule__ElkSingleEdgeSection__UnorderedGroup_1_0__3
                    {
                    pushFollow(FOLLOW_2);
                    rule__ElkSingleEdgeSection__UnorderedGroup_1_0__3();

                    state._fsp--;


                    }
                    break;

            }


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ElkSingleEdgeSection__UnorderedGroup_1_0__2"


    // $ANTLR start "rule__ElkSingleEdgeSection__UnorderedGroup_1_0__3"
    // InternalElkGraph.g:4898:1: rule__ElkSingleEdgeSection__UnorderedGroup_1_0__3 : rule__ElkSingleEdgeSection__UnorderedGroup_1_0__Impl ;
    public final void rule__ElkSingleEdgeSection__UnorderedGroup_1_0__3() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:4902:1: ( rule__ElkSingleEdgeSection__UnorderedGroup_1_0__Impl )
            // InternalElkGraph.g:4903:2: rule__ElkSingleEdgeSection__UnorderedGroup_1_0__Impl
            {
            pushFollow(FOLLOW_2);
            rule__ElkSingleEdgeSection__UnorderedGroup_1_0__Impl();

            state._fsp--;


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ElkSingleEdgeSection__UnorderedGroup_1_0__3"


    // $ANTLR start "rule__ElkEdgeSection__UnorderedGroup_4_0"
    // InternalElkGraph.g:4910:1: rule__ElkEdgeSection__UnorderedGroup_4_0 : ( rule__ElkEdgeSection__UnorderedGroup_4_0__0 )? ;
    public final void rule__ElkEdgeSection__UnorderedGroup_4_0() throws RecognitionException {

        		int stackSize = keepStackSize();
        		getUnorderedGroupHelper().enter(grammarAccess.getElkEdgeSectionAccess().getUnorderedGroup_4_0());
        	
        try {
            // InternalElkGraph.g:4915:1: ( ( rule__ElkEdgeSection__UnorderedGroup_4_0__0 )? )
            // InternalElkGraph.g:4916:2: ( rule__ElkEdgeSection__UnorderedGroup_4_0__0 )?
            {
            // InternalElkGraph.g:4916:2: ( rule__ElkEdgeSection__UnorderedGroup_4_0__0 )?
            int alt51=2;
            int LA51_0 = input.LA(1);

            if ( LA51_0 == 31 && getUnorderedGroupHelper().canSelect(grammarAccess.getElkEdgeSectionAccess().getUnorderedGroup_4_0(), 0) ) {
                alt51=1;
            }
            else if ( LA51_0 == 32 && getUnorderedGroupHelper().canSelect(grammarAccess.getElkEdgeSectionAccess().getUnorderedGroup_4_0(), 1) ) {
                alt51=1;
            }
            else if ( LA51_0 == 33 && getUnorderedGroupHelper().canSelect(grammarAccess.getElkEdgeSectionAccess().getUnorderedGroup_4_0(), 2) ) {
                alt51=1;
            }
            else if ( LA51_0 == 34 && getUnorderedGroupHelper().canSelect(grammarAccess.getElkEdgeSectionAccess().getUnorderedGroup_4_0(), 3) ) {
                alt51=1;
            }
            switch (alt51) {
                case 1 :
                    // InternalElkGraph.g:4916:2: rule__ElkEdgeSection__UnorderedGroup_4_0__0
                    {
                    pushFollow(FOLLOW_2);
                    rule__ElkEdgeSection__UnorderedGroup_4_0__0();

                    state._fsp--;


                    }
                    break;

            }


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	getUnorderedGroupHelper().leave(grammarAccess.getElkEdgeSectionAccess().getUnorderedGroup_4_0());
            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ElkEdgeSection__UnorderedGroup_4_0"


    // $ANTLR start "rule__ElkEdgeSection__UnorderedGroup_4_0__Impl"
    // InternalElkGraph.g:4924:1: rule__ElkEdgeSection__UnorderedGroup_4_0__Impl : ( ({...}? => ( ( ( rule__ElkEdgeSection__Group_4_0_0__0 ) ) ) ) | ({...}? => ( ( ( rule__ElkEdgeSection__Group_4_0_1__0 ) ) ) ) | ({...}? => ( ( ( rule__ElkEdgeSection__Group_4_0_2__0 ) ) ) ) | ({...}? => ( ( ( rule__ElkEdgeSection__Group_4_0_3__0 ) ) ) ) ) ;
    public final void rule__ElkEdgeSection__UnorderedGroup_4_0__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        		boolean selected = false;
        	
        try {
            // InternalElkGraph.g:4929:1: ( ( ({...}? => ( ( ( rule__ElkEdgeSection__Group_4_0_0__0 ) ) ) ) | ({...}? => ( ( ( rule__ElkEdgeSection__Group_4_0_1__0 ) ) ) ) | ({...}? => ( ( ( rule__ElkEdgeSection__Group_4_0_2__0 ) ) ) ) | ({...}? => ( ( ( rule__ElkEdgeSection__Group_4_0_3__0 ) ) ) ) ) )
            // InternalElkGraph.g:4930:3: ( ({...}? => ( ( ( rule__ElkEdgeSection__Group_4_0_0__0 ) ) ) ) | ({...}? => ( ( ( rule__ElkEdgeSection__Group_4_0_1__0 ) ) ) ) | ({...}? => ( ( ( rule__ElkEdgeSection__Group_4_0_2__0 ) ) ) ) | ({...}? => ( ( ( rule__ElkEdgeSection__Group_4_0_3__0 ) ) ) ) )
            {
            // InternalElkGraph.g:4930:3: ( ({...}? => ( ( ( rule__ElkEdgeSection__Group_4_0_0__0 ) ) ) ) | ({...}? => ( ( ( rule__ElkEdgeSection__Group_4_0_1__0 ) ) ) ) | ({...}? => ( ( ( rule__ElkEdgeSection__Group_4_0_2__0 ) ) ) ) | ({...}? => ( ( ( rule__ElkEdgeSection__Group_4_0_3__0 ) ) ) ) )
            int alt52=4;
            int LA52_0 = input.LA(1);

            if ( LA52_0 == 31 && getUnorderedGroupHelper().canSelect(grammarAccess.getElkEdgeSectionAccess().getUnorderedGroup_4_0(), 0) ) {
                alt52=1;
            }
            else if ( LA52_0 == 32 && getUnorderedGroupHelper().canSelect(grammarAccess.getElkEdgeSectionAccess().getUnorderedGroup_4_0(), 1) ) {
                alt52=2;
            }
            else if ( LA52_0 == 33 && getUnorderedGroupHelper().canSelect(grammarAccess.getElkEdgeSectionAccess().getUnorderedGroup_4_0(), 2) ) {
                alt52=3;
            }
            else if ( LA52_0 == 34 && getUnorderedGroupHelper().canSelect(grammarAccess.getElkEdgeSectionAccess().getUnorderedGroup_4_0(), 3) ) {
                alt52=4;
            }
            else {
                NoViableAltException nvae =
                    new NoViableAltException("", 52, 0, input);

                throw nvae;
            }
            switch (alt52) {
                case 1 :
                    // InternalElkGraph.g:4931:3: ({...}? => ( ( ( rule__ElkEdgeSection__Group_4_0_0__0 ) ) ) )
                    {
                    // InternalElkGraph.g:4931:3: ({...}? => ( ( ( rule__ElkEdgeSection__Group_4_0_0__0 ) ) ) )
                    // InternalElkGraph.g:4932:4: {...}? => ( ( ( rule__ElkEdgeSection__Group_4_0_0__0 ) ) )
                    {
                    if ( ! getUnorderedGroupHelper().canSelect(grammarAccess.getElkEdgeSectionAccess().getUnorderedGroup_4_0(), 0) ) {
                        throw new FailedPredicateException(input, "rule__ElkEdgeSection__UnorderedGroup_4_0__Impl", "getUnorderedGroupHelper().canSelect(grammarAccess.getElkEdgeSectionAccess().getUnorderedGroup_4_0(), 0)");
                    }
                    // InternalElkGraph.g:4932:112: ( ( ( rule__ElkEdgeSection__Group_4_0_0__0 ) ) )
                    // InternalElkGraph.g:4933:5: ( ( rule__ElkEdgeSection__Group_4_0_0__0 ) )
                    {

                    					getUnorderedGroupHelper().select(grammarAccess.getElkEdgeSectionAccess().getUnorderedGroup_4_0(), 0);
                    				

                    					selected = true;
                    				
                    // InternalElkGraph.g:4939:5: ( ( rule__ElkEdgeSection__Group_4_0_0__0 ) )
                    // InternalElkGraph.g:4940:6: ( rule__ElkEdgeSection__Group_4_0_0__0 )
                    {
                     before(grammarAccess.getElkEdgeSectionAccess().getGroup_4_0_0()); 
                    // InternalElkGraph.g:4941:6: ( rule__ElkEdgeSection__Group_4_0_0__0 )
                    // InternalElkGraph.g:4941:7: rule__ElkEdgeSection__Group_4_0_0__0
                    {
                    pushFollow(FOLLOW_2);
                    rule__ElkEdgeSection__Group_4_0_0__0();

                    state._fsp--;


                    }

                     after(grammarAccess.getElkEdgeSectionAccess().getGroup_4_0_0()); 

                    }


                    }


                    }


                    }
                    break;
                case 2 :
                    // InternalElkGraph.g:4946:3: ({...}? => ( ( ( rule__ElkEdgeSection__Group_4_0_1__0 ) ) ) )
                    {
                    // InternalElkGraph.g:4946:3: ({...}? => ( ( ( rule__ElkEdgeSection__Group_4_0_1__0 ) ) ) )
                    // InternalElkGraph.g:4947:4: {...}? => ( ( ( rule__ElkEdgeSection__Group_4_0_1__0 ) ) )
                    {
                    if ( ! getUnorderedGroupHelper().canSelect(grammarAccess.getElkEdgeSectionAccess().getUnorderedGroup_4_0(), 1) ) {
                        throw new FailedPredicateException(input, "rule__ElkEdgeSection__UnorderedGroup_4_0__Impl", "getUnorderedGroupHelper().canSelect(grammarAccess.getElkEdgeSectionAccess().getUnorderedGroup_4_0(), 1)");
                    }
                    // InternalElkGraph.g:4947:112: ( ( ( rule__ElkEdgeSection__Group_4_0_1__0 ) ) )
                    // InternalElkGraph.g:4948:5: ( ( rule__ElkEdgeSection__Group_4_0_1__0 ) )
                    {

                    					getUnorderedGroupHelper().select(grammarAccess.getElkEdgeSectionAccess().getUnorderedGroup_4_0(), 1);
                    				

                    					selected = true;
                    				
                    // InternalElkGraph.g:4954:5: ( ( rule__ElkEdgeSection__Group_4_0_1__0 ) )
                    // InternalElkGraph.g:4955:6: ( rule__ElkEdgeSection__Group_4_0_1__0 )
                    {
                     before(grammarAccess.getElkEdgeSectionAccess().getGroup_4_0_1()); 
                    // InternalElkGraph.g:4956:6: ( rule__ElkEdgeSection__Group_4_0_1__0 )
                    // InternalElkGraph.g:4956:7: rule__ElkEdgeSection__Group_4_0_1__0
                    {
                    pushFollow(FOLLOW_2);
                    rule__ElkEdgeSection__Group_4_0_1__0();

                    state._fsp--;


                    }

                     after(grammarAccess.getElkEdgeSectionAccess().getGroup_4_0_1()); 

                    }


                    }


                    }


                    }
                    break;
                case 3 :
                    // InternalElkGraph.g:4961:3: ({...}? => ( ( ( rule__ElkEdgeSection__Group_4_0_2__0 ) ) ) )
                    {
                    // InternalElkGraph.g:4961:3: ({...}? => ( ( ( rule__ElkEdgeSection__Group_4_0_2__0 ) ) ) )
                    // InternalElkGraph.g:4962:4: {...}? => ( ( ( rule__ElkEdgeSection__Group_4_0_2__0 ) ) )
                    {
                    if ( ! getUnorderedGroupHelper().canSelect(grammarAccess.getElkEdgeSectionAccess().getUnorderedGroup_4_0(), 2) ) {
                        throw new FailedPredicateException(input, "rule__ElkEdgeSection__UnorderedGroup_4_0__Impl", "getUnorderedGroupHelper().canSelect(grammarAccess.getElkEdgeSectionAccess().getUnorderedGroup_4_0(), 2)");
                    }
                    // InternalElkGraph.g:4962:112: ( ( ( rule__ElkEdgeSection__Group_4_0_2__0 ) ) )
                    // InternalElkGraph.g:4963:5: ( ( rule__ElkEdgeSection__Group_4_0_2__0 ) )
                    {

                    					getUnorderedGroupHelper().select(grammarAccess.getElkEdgeSectionAccess().getUnorderedGroup_4_0(), 2);
                    				

                    					selected = true;
                    				
                    // InternalElkGraph.g:4969:5: ( ( rule__ElkEdgeSection__Group_4_0_2__0 ) )
                    // InternalElkGraph.g:4970:6: ( rule__ElkEdgeSection__Group_4_0_2__0 )
                    {
                     before(grammarAccess.getElkEdgeSectionAccess().getGroup_4_0_2()); 
                    // InternalElkGraph.g:4971:6: ( rule__ElkEdgeSection__Group_4_0_2__0 )
                    // InternalElkGraph.g:4971:7: rule__ElkEdgeSection__Group_4_0_2__0
                    {
                    pushFollow(FOLLOW_2);
                    rule__ElkEdgeSection__Group_4_0_2__0();

                    state._fsp--;


                    }

                     after(grammarAccess.getElkEdgeSectionAccess().getGroup_4_0_2()); 

                    }


                    }


                    }


                    }
                    break;
                case 4 :
                    // InternalElkGraph.g:4976:3: ({...}? => ( ( ( rule__ElkEdgeSection__Group_4_0_3__0 ) ) ) )
                    {
                    // InternalElkGraph.g:4976:3: ({...}? => ( ( ( rule__ElkEdgeSection__Group_4_0_3__0 ) ) ) )
                    // InternalElkGraph.g:4977:4: {...}? => ( ( ( rule__ElkEdgeSection__Group_4_0_3__0 ) ) )
                    {
                    if ( ! getUnorderedGroupHelper().canSelect(grammarAccess.getElkEdgeSectionAccess().getUnorderedGroup_4_0(), 3) ) {
                        throw new FailedPredicateException(input, "rule__ElkEdgeSection__UnorderedGroup_4_0__Impl", "getUnorderedGroupHelper().canSelect(grammarAccess.getElkEdgeSectionAccess().getUnorderedGroup_4_0(), 3)");
                    }
                    // InternalElkGraph.g:4977:112: ( ( ( rule__ElkEdgeSection__Group_4_0_3__0 ) ) )
                    // InternalElkGraph.g:4978:5: ( ( rule__ElkEdgeSection__Group_4_0_3__0 ) )
                    {

                    					getUnorderedGroupHelper().select(grammarAccess.getElkEdgeSectionAccess().getUnorderedGroup_4_0(), 3);
                    				

                    					selected = true;
                    				
                    // InternalElkGraph.g:4984:5: ( ( rule__ElkEdgeSection__Group_4_0_3__0 ) )
                    // InternalElkGraph.g:4985:6: ( rule__ElkEdgeSection__Group_4_0_3__0 )
                    {
                     before(grammarAccess.getElkEdgeSectionAccess().getGroup_4_0_3()); 
                    // InternalElkGraph.g:4986:6: ( rule__ElkEdgeSection__Group_4_0_3__0 )
                    // InternalElkGraph.g:4986:7: rule__ElkEdgeSection__Group_4_0_3__0
                    {
                    pushFollow(FOLLOW_2);
                    rule__ElkEdgeSection__Group_4_0_3__0();

                    state._fsp--;


                    }

                     after(grammarAccess.getElkEdgeSectionAccess().getGroup_4_0_3()); 

                    }


                    }


                    }


                    }
                    break;

            }


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	if (selected)
            		getUnorderedGroupHelper().returnFromSelection(grammarAccess.getElkEdgeSectionAccess().getUnorderedGroup_4_0());
            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ElkEdgeSection__UnorderedGroup_4_0__Impl"


    // $ANTLR start "rule__ElkEdgeSection__UnorderedGroup_4_0__0"
    // InternalElkGraph.g:4999:1: rule__ElkEdgeSection__UnorderedGroup_4_0__0 : rule__ElkEdgeSection__UnorderedGroup_4_0__Impl ( rule__ElkEdgeSection__UnorderedGroup_4_0__1 )? ;
    public final void rule__ElkEdgeSection__UnorderedGroup_4_0__0() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:5003:1: ( rule__ElkEdgeSection__UnorderedGroup_4_0__Impl ( rule__ElkEdgeSection__UnorderedGroup_4_0__1 )? )
            // InternalElkGraph.g:5004:2: rule__ElkEdgeSection__UnorderedGroup_4_0__Impl ( rule__ElkEdgeSection__UnorderedGroup_4_0__1 )?
            {
            pushFollow(FOLLOW_32);
            rule__ElkEdgeSection__UnorderedGroup_4_0__Impl();

            state._fsp--;

            // InternalElkGraph.g:5005:2: ( rule__ElkEdgeSection__UnorderedGroup_4_0__1 )?
            int alt53=2;
            int LA53_0 = input.LA(1);

            if ( LA53_0 == 31 && getUnorderedGroupHelper().canSelect(grammarAccess.getElkEdgeSectionAccess().getUnorderedGroup_4_0(), 0) ) {
                alt53=1;
            }
            else if ( LA53_0 == 32 && getUnorderedGroupHelper().canSelect(grammarAccess.getElkEdgeSectionAccess().getUnorderedGroup_4_0(), 1) ) {
                alt53=1;
            }
            else if ( LA53_0 == 33 && getUnorderedGroupHelper().canSelect(grammarAccess.getElkEdgeSectionAccess().getUnorderedGroup_4_0(), 2) ) {
                alt53=1;
            }
            else if ( LA53_0 == 34 && getUnorderedGroupHelper().canSelect(grammarAccess.getElkEdgeSectionAccess().getUnorderedGroup_4_0(), 3) ) {
                alt53=1;
            }
            switch (alt53) {
                case 1 :
                    // InternalElkGraph.g:5005:2: rule__ElkEdgeSection__UnorderedGroup_4_0__1
                    {
                    pushFollow(FOLLOW_2);
                    rule__ElkEdgeSection__UnorderedGroup_4_0__1();

                    state._fsp--;


                    }
                    break;

            }


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ElkEdgeSection__UnorderedGroup_4_0__0"


    // $ANTLR start "rule__ElkEdgeSection__UnorderedGroup_4_0__1"
    // InternalElkGraph.g:5011:1: rule__ElkEdgeSection__UnorderedGroup_4_0__1 : rule__ElkEdgeSection__UnorderedGroup_4_0__Impl ( rule__ElkEdgeSection__UnorderedGroup_4_0__2 )? ;
    public final void rule__ElkEdgeSection__UnorderedGroup_4_0__1() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:5015:1: ( rule__ElkEdgeSection__UnorderedGroup_4_0__Impl ( rule__ElkEdgeSection__UnorderedGroup_4_0__2 )? )
            // InternalElkGraph.g:5016:2: rule__ElkEdgeSection__UnorderedGroup_4_0__Impl ( rule__ElkEdgeSection__UnorderedGroup_4_0__2 )?
            {
            pushFollow(FOLLOW_32);
            rule__ElkEdgeSection__UnorderedGroup_4_0__Impl();

            state._fsp--;

            // InternalElkGraph.g:5017:2: ( rule__ElkEdgeSection__UnorderedGroup_4_0__2 )?
            int alt54=2;
            int LA54_0 = input.LA(1);

            if ( LA54_0 == 31 && getUnorderedGroupHelper().canSelect(grammarAccess.getElkEdgeSectionAccess().getUnorderedGroup_4_0(), 0) ) {
                alt54=1;
            }
            else if ( LA54_0 == 32 && getUnorderedGroupHelper().canSelect(grammarAccess.getElkEdgeSectionAccess().getUnorderedGroup_4_0(), 1) ) {
                alt54=1;
            }
            else if ( LA54_0 == 33 && getUnorderedGroupHelper().canSelect(grammarAccess.getElkEdgeSectionAccess().getUnorderedGroup_4_0(), 2) ) {
                alt54=1;
            }
            else if ( LA54_0 == 34 && getUnorderedGroupHelper().canSelect(grammarAccess.getElkEdgeSectionAccess().getUnorderedGroup_4_0(), 3) ) {
                alt54=1;
            }
            switch (alt54) {
                case 1 :
                    // InternalElkGraph.g:5017:2: rule__ElkEdgeSection__UnorderedGroup_4_0__2
                    {
                    pushFollow(FOLLOW_2);
                    rule__ElkEdgeSection__UnorderedGroup_4_0__2();

                    state._fsp--;


                    }
                    break;

            }


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ElkEdgeSection__UnorderedGroup_4_0__1"


    // $ANTLR start "rule__ElkEdgeSection__UnorderedGroup_4_0__2"
    // InternalElkGraph.g:5023:1: rule__ElkEdgeSection__UnorderedGroup_4_0__2 : rule__ElkEdgeSection__UnorderedGroup_4_0__Impl ( rule__ElkEdgeSection__UnorderedGroup_4_0__3 )? ;
    public final void rule__ElkEdgeSection__UnorderedGroup_4_0__2() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:5027:1: ( rule__ElkEdgeSection__UnorderedGroup_4_0__Impl ( rule__ElkEdgeSection__UnorderedGroup_4_0__3 )? )
            // InternalElkGraph.g:5028:2: rule__ElkEdgeSection__UnorderedGroup_4_0__Impl ( rule__ElkEdgeSection__UnorderedGroup_4_0__3 )?
            {
            pushFollow(FOLLOW_32);
            rule__ElkEdgeSection__UnorderedGroup_4_0__Impl();

            state._fsp--;

            // InternalElkGraph.g:5029:2: ( rule__ElkEdgeSection__UnorderedGroup_4_0__3 )?
            int alt55=2;
            int LA55_0 = input.LA(1);

            if ( LA55_0 == 31 && getUnorderedGroupHelper().canSelect(grammarAccess.getElkEdgeSectionAccess().getUnorderedGroup_4_0(), 0) ) {
                alt55=1;
            }
            else if ( LA55_0 == 32 && getUnorderedGroupHelper().canSelect(grammarAccess.getElkEdgeSectionAccess().getUnorderedGroup_4_0(), 1) ) {
                alt55=1;
            }
            else if ( LA55_0 == 33 && getUnorderedGroupHelper().canSelect(grammarAccess.getElkEdgeSectionAccess().getUnorderedGroup_4_0(), 2) ) {
                alt55=1;
            }
            else if ( LA55_0 == 34 && getUnorderedGroupHelper().canSelect(grammarAccess.getElkEdgeSectionAccess().getUnorderedGroup_4_0(), 3) ) {
                alt55=1;
            }
            switch (alt55) {
                case 1 :
                    // InternalElkGraph.g:5029:2: rule__ElkEdgeSection__UnorderedGroup_4_0__3
                    {
                    pushFollow(FOLLOW_2);
                    rule__ElkEdgeSection__UnorderedGroup_4_0__3();

                    state._fsp--;


                    }
                    break;

            }


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ElkEdgeSection__UnorderedGroup_4_0__2"


    // $ANTLR start "rule__ElkEdgeSection__UnorderedGroup_4_0__3"
    // InternalElkGraph.g:5035:1: rule__ElkEdgeSection__UnorderedGroup_4_0__3 : rule__ElkEdgeSection__UnorderedGroup_4_0__Impl ;
    public final void rule__ElkEdgeSection__UnorderedGroup_4_0__3() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:5039:1: ( rule__ElkEdgeSection__UnorderedGroup_4_0__Impl )
            // InternalElkGraph.g:5040:2: rule__ElkEdgeSection__UnorderedGroup_4_0__Impl
            {
            pushFollow(FOLLOW_2);
            rule__ElkEdgeSection__UnorderedGroup_4_0__Impl();

            state._fsp--;


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ElkEdgeSection__UnorderedGroup_4_0__3"


    // $ANTLR start "rule__RootNode__IdentifierAssignment_1_1"
    // InternalElkGraph.g:5047:1: rule__RootNode__IdentifierAssignment_1_1 : ( RULE_ID ) ;
    public final void rule__RootNode__IdentifierAssignment_1_1() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:5051:1: ( ( RULE_ID ) )
            // InternalElkGraph.g:5052:2: ( RULE_ID )
            {
            // InternalElkGraph.g:5052:2: ( RULE_ID )
            // InternalElkGraph.g:5053:3: RULE_ID
            {
             before(grammarAccess.getRootNodeAccess().getIdentifierIDTerminalRuleCall_1_1_0()); 
            match(input,RULE_ID,FOLLOW_2); 
             after(grammarAccess.getRootNodeAccess().getIdentifierIDTerminalRuleCall_1_1_0()); 

            }


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__RootNode__IdentifierAssignment_1_1"


    // $ANTLR start "rule__RootNode__PropertiesAssignment_3"
    // InternalElkGraph.g:5062:1: rule__RootNode__PropertiesAssignment_3 : ( ruleProperty ) ;
    public final void rule__RootNode__PropertiesAssignment_3() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:5066:1: ( ( ruleProperty ) )
            // InternalElkGraph.g:5067:2: ( ruleProperty )
            {
            // InternalElkGraph.g:5067:2: ( ruleProperty )
            // InternalElkGraph.g:5068:3: ruleProperty
            {
             before(grammarAccess.getRootNodeAccess().getPropertiesPropertyParserRuleCall_3_0()); 
            pushFollow(FOLLOW_2);
            ruleProperty();

            state._fsp--;

             after(grammarAccess.getRootNodeAccess().getPropertiesPropertyParserRuleCall_3_0()); 

            }


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__RootNode__PropertiesAssignment_3"


    // $ANTLR start "rule__RootNode__LabelsAssignment_4_0"
    // InternalElkGraph.g:5077:1: rule__RootNode__LabelsAssignment_4_0 : ( ruleElkLabel ) ;
    public final void rule__RootNode__LabelsAssignment_4_0() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:5081:1: ( ( ruleElkLabel ) )
            // InternalElkGraph.g:5082:2: ( ruleElkLabel )
            {
            // InternalElkGraph.g:5082:2: ( ruleElkLabel )
            // InternalElkGraph.g:5083:3: ruleElkLabel
            {
             before(grammarAccess.getRootNodeAccess().getLabelsElkLabelParserRuleCall_4_0_0()); 
            pushFollow(FOLLOW_2);
            ruleElkLabel();

            state._fsp--;

             after(grammarAccess.getRootNodeAccess().getLabelsElkLabelParserRuleCall_4_0_0()); 

            }


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__RootNode__LabelsAssignment_4_0"


    // $ANTLR start "rule__RootNode__PortsAssignment_4_1"
    // InternalElkGraph.g:5092:1: rule__RootNode__PortsAssignment_4_1 : ( ruleElkPort ) ;
    public final void rule__RootNode__PortsAssignment_4_1() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:5096:1: ( ( ruleElkPort ) )
            // InternalElkGraph.g:5097:2: ( ruleElkPort )
            {
            // InternalElkGraph.g:5097:2: ( ruleElkPort )
            // InternalElkGraph.g:5098:3: ruleElkPort
            {
             before(grammarAccess.getRootNodeAccess().getPortsElkPortParserRuleCall_4_1_0()); 
            pushFollow(FOLLOW_2);
            ruleElkPort();

            state._fsp--;

             after(grammarAccess.getRootNodeAccess().getPortsElkPortParserRuleCall_4_1_0()); 

            }


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__RootNode__PortsAssignment_4_1"


    // $ANTLR start "rule__RootNode__ChildrenAssignment_4_2"
    // InternalElkGraph.g:5107:1: rule__RootNode__ChildrenAssignment_4_2 : ( ruleElkNode ) ;
    public final void rule__RootNode__ChildrenAssignment_4_2() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:5111:1: ( ( ruleElkNode ) )
            // InternalElkGraph.g:5112:2: ( ruleElkNode )
            {
            // InternalElkGraph.g:5112:2: ( ruleElkNode )
            // InternalElkGraph.g:5113:3: ruleElkNode
            {
             before(grammarAccess.getRootNodeAccess().getChildrenElkNodeParserRuleCall_4_2_0()); 
            pushFollow(FOLLOW_2);
            ruleElkNode();

            state._fsp--;

             after(grammarAccess.getRootNodeAccess().getChildrenElkNodeParserRuleCall_4_2_0()); 

            }


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__RootNode__ChildrenAssignment_4_2"


    // $ANTLR start "rule__RootNode__ContainedEdgesAssignment_4_3"
    // InternalElkGraph.g:5122:1: rule__RootNode__ContainedEdgesAssignment_4_3 : ( ruleElkEdge ) ;
    public final void rule__RootNode__ContainedEdgesAssignment_4_3() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:5126:1: ( ( ruleElkEdge ) )
            // InternalElkGraph.g:5127:2: ( ruleElkEdge )
            {
            // InternalElkGraph.g:5127:2: ( ruleElkEdge )
            // InternalElkGraph.g:5128:3: ruleElkEdge
            {
             before(grammarAccess.getRootNodeAccess().getContainedEdgesElkEdgeParserRuleCall_4_3_0()); 
            pushFollow(FOLLOW_2);
            ruleElkEdge();

            state._fsp--;

             after(grammarAccess.getRootNodeAccess().getContainedEdgesElkEdgeParserRuleCall_4_3_0()); 

            }


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__RootNode__ContainedEdgesAssignment_4_3"


    // $ANTLR start "rule__ElkNode__IdentifierAssignment_1"
    // InternalElkGraph.g:5137:1: rule__ElkNode__IdentifierAssignment_1 : ( RULE_ID ) ;
    public final void rule__ElkNode__IdentifierAssignment_1() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:5141:1: ( ( RULE_ID ) )
            // InternalElkGraph.g:5142:2: ( RULE_ID )
            {
            // InternalElkGraph.g:5142:2: ( RULE_ID )
            // InternalElkGraph.g:5143:3: RULE_ID
            {
             before(grammarAccess.getElkNodeAccess().getIdentifierIDTerminalRuleCall_1_0()); 
            match(input,RULE_ID,FOLLOW_2); 
             after(grammarAccess.getElkNodeAccess().getIdentifierIDTerminalRuleCall_1_0()); 

            }


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ElkNode__IdentifierAssignment_1"


    // $ANTLR start "rule__ElkNode__PropertiesAssignment_2_2"
    // InternalElkGraph.g:5152:1: rule__ElkNode__PropertiesAssignment_2_2 : ( ruleProperty ) ;
    public final void rule__ElkNode__PropertiesAssignment_2_2() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:5156:1: ( ( ruleProperty ) )
            // InternalElkGraph.g:5157:2: ( ruleProperty )
            {
            // InternalElkGraph.g:5157:2: ( ruleProperty )
            // InternalElkGraph.g:5158:3: ruleProperty
            {
             before(grammarAccess.getElkNodeAccess().getPropertiesPropertyParserRuleCall_2_2_0()); 
            pushFollow(FOLLOW_2);
            ruleProperty();

            state._fsp--;

             after(grammarAccess.getElkNodeAccess().getPropertiesPropertyParserRuleCall_2_2_0()); 

            }


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ElkNode__PropertiesAssignment_2_2"


    // $ANTLR start "rule__ElkNode__LabelsAssignment_2_3_0"
    // InternalElkGraph.g:5167:1: rule__ElkNode__LabelsAssignment_2_3_0 : ( ruleElkLabel ) ;
    public final void rule__ElkNode__LabelsAssignment_2_3_0() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:5171:1: ( ( ruleElkLabel ) )
            // InternalElkGraph.g:5172:2: ( ruleElkLabel )
            {
            // InternalElkGraph.g:5172:2: ( ruleElkLabel )
            // InternalElkGraph.g:5173:3: ruleElkLabel
            {
             before(grammarAccess.getElkNodeAccess().getLabelsElkLabelParserRuleCall_2_3_0_0()); 
            pushFollow(FOLLOW_2);
            ruleElkLabel();

            state._fsp--;

             after(grammarAccess.getElkNodeAccess().getLabelsElkLabelParserRuleCall_2_3_0_0()); 

            }


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ElkNode__LabelsAssignment_2_3_0"


    // $ANTLR start "rule__ElkNode__PortsAssignment_2_3_1"
    // InternalElkGraph.g:5182:1: rule__ElkNode__PortsAssignment_2_3_1 : ( ruleElkPort ) ;
    public final void rule__ElkNode__PortsAssignment_2_3_1() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:5186:1: ( ( ruleElkPort ) )
            // InternalElkGraph.g:5187:2: ( ruleElkPort )
            {
            // InternalElkGraph.g:5187:2: ( ruleElkPort )
            // InternalElkGraph.g:5188:3: ruleElkPort
            {
             before(grammarAccess.getElkNodeAccess().getPortsElkPortParserRuleCall_2_3_1_0()); 
            pushFollow(FOLLOW_2);
            ruleElkPort();

            state._fsp--;

             after(grammarAccess.getElkNodeAccess().getPortsElkPortParserRuleCall_2_3_1_0()); 

            }


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ElkNode__PortsAssignment_2_3_1"


    // $ANTLR start "rule__ElkNode__ChildrenAssignment_2_3_2"
    // InternalElkGraph.g:5197:1: rule__ElkNode__ChildrenAssignment_2_3_2 : ( ruleElkNode ) ;
    public final void rule__ElkNode__ChildrenAssignment_2_3_2() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:5201:1: ( ( ruleElkNode ) )
            // InternalElkGraph.g:5202:2: ( ruleElkNode )
            {
            // InternalElkGraph.g:5202:2: ( ruleElkNode )
            // InternalElkGraph.g:5203:3: ruleElkNode
            {
             before(grammarAccess.getElkNodeAccess().getChildrenElkNodeParserRuleCall_2_3_2_0()); 
            pushFollow(FOLLOW_2);
            ruleElkNode();

            state._fsp--;

             after(grammarAccess.getElkNodeAccess().getChildrenElkNodeParserRuleCall_2_3_2_0()); 

            }


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ElkNode__ChildrenAssignment_2_3_2"


    // $ANTLR start "rule__ElkNode__ContainedEdgesAssignment_2_3_3"
    // InternalElkGraph.g:5212:1: rule__ElkNode__ContainedEdgesAssignment_2_3_3 : ( ruleElkEdge ) ;
    public final void rule__ElkNode__ContainedEdgesAssignment_2_3_3() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:5216:1: ( ( ruleElkEdge ) )
            // InternalElkGraph.g:5217:2: ( ruleElkEdge )
            {
            // InternalElkGraph.g:5217:2: ( ruleElkEdge )
            // InternalElkGraph.g:5218:3: ruleElkEdge
            {
             before(grammarAccess.getElkNodeAccess().getContainedEdgesElkEdgeParserRuleCall_2_3_3_0()); 
            pushFollow(FOLLOW_2);
            ruleElkEdge();

            state._fsp--;

             after(grammarAccess.getElkNodeAccess().getContainedEdgesElkEdgeParserRuleCall_2_3_3_0()); 

            }


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ElkNode__ContainedEdgesAssignment_2_3_3"


    // $ANTLR start "rule__ElkLabel__IdentifierAssignment_1_0"
    // InternalElkGraph.g:5227:1: rule__ElkLabel__IdentifierAssignment_1_0 : ( RULE_ID ) ;
    public final void rule__ElkLabel__IdentifierAssignment_1_0() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:5231:1: ( ( RULE_ID ) )
            // InternalElkGraph.g:5232:2: ( RULE_ID )
            {
            // InternalElkGraph.g:5232:2: ( RULE_ID )
            // InternalElkGraph.g:5233:3: RULE_ID
            {
             before(grammarAccess.getElkLabelAccess().getIdentifierIDTerminalRuleCall_1_0_0()); 
            match(input,RULE_ID,FOLLOW_2); 
             after(grammarAccess.getElkLabelAccess().getIdentifierIDTerminalRuleCall_1_0_0()); 

            }


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ElkLabel__IdentifierAssignment_1_0"


    // $ANTLR start "rule__ElkLabel__TextAssignment_2"
    // InternalElkGraph.g:5242:1: rule__ElkLabel__TextAssignment_2 : ( RULE_STRING ) ;
    public final void rule__ElkLabel__TextAssignment_2() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:5246:1: ( ( RULE_STRING ) )
            // InternalElkGraph.g:5247:2: ( RULE_STRING )
            {
            // InternalElkGraph.g:5247:2: ( RULE_STRING )
            // InternalElkGraph.g:5248:3: RULE_STRING
            {
             before(grammarAccess.getElkLabelAccess().getTextSTRINGTerminalRuleCall_2_0()); 
            match(input,RULE_STRING,FOLLOW_2); 
             after(grammarAccess.getElkLabelAccess().getTextSTRINGTerminalRuleCall_2_0()); 

            }


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ElkLabel__TextAssignment_2"


    // $ANTLR start "rule__ElkLabel__PropertiesAssignment_3_2"
    // InternalElkGraph.g:5257:1: rule__ElkLabel__PropertiesAssignment_3_2 : ( ruleProperty ) ;
    public final void rule__ElkLabel__PropertiesAssignment_3_2() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:5261:1: ( ( ruleProperty ) )
            // InternalElkGraph.g:5262:2: ( ruleProperty )
            {
            // InternalElkGraph.g:5262:2: ( ruleProperty )
            // InternalElkGraph.g:5263:3: ruleProperty
            {
             before(grammarAccess.getElkLabelAccess().getPropertiesPropertyParserRuleCall_3_2_0()); 
            pushFollow(FOLLOW_2);
            ruleProperty();

            state._fsp--;

             after(grammarAccess.getElkLabelAccess().getPropertiesPropertyParserRuleCall_3_2_0()); 

            }


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ElkLabel__PropertiesAssignment_3_2"


    // $ANTLR start "rule__ElkLabel__LabelsAssignment_3_3"
    // InternalElkGraph.g:5272:1: rule__ElkLabel__LabelsAssignment_3_3 : ( ruleElkLabel ) ;
    public final void rule__ElkLabel__LabelsAssignment_3_3() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:5276:1: ( ( ruleElkLabel ) )
            // InternalElkGraph.g:5277:2: ( ruleElkLabel )
            {
            // InternalElkGraph.g:5277:2: ( ruleElkLabel )
            // InternalElkGraph.g:5278:3: ruleElkLabel
            {
             before(grammarAccess.getElkLabelAccess().getLabelsElkLabelParserRuleCall_3_3_0()); 
            pushFollow(FOLLOW_2);
            ruleElkLabel();

            state._fsp--;

             after(grammarAccess.getElkLabelAccess().getLabelsElkLabelParserRuleCall_3_3_0()); 

            }


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ElkLabel__LabelsAssignment_3_3"


    // $ANTLR start "rule__ElkPort__IdentifierAssignment_1"
    // InternalElkGraph.g:5287:1: rule__ElkPort__IdentifierAssignment_1 : ( RULE_ID ) ;
    public final void rule__ElkPort__IdentifierAssignment_1() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:5291:1: ( ( RULE_ID ) )
            // InternalElkGraph.g:5292:2: ( RULE_ID )
            {
            // InternalElkGraph.g:5292:2: ( RULE_ID )
            // InternalElkGraph.g:5293:3: RULE_ID
            {
             before(grammarAccess.getElkPortAccess().getIdentifierIDTerminalRuleCall_1_0()); 
            match(input,RULE_ID,FOLLOW_2); 
             after(grammarAccess.getElkPortAccess().getIdentifierIDTerminalRuleCall_1_0()); 

            }


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ElkPort__IdentifierAssignment_1"


    // $ANTLR start "rule__ElkPort__PropertiesAssignment_2_2"
    // InternalElkGraph.g:5302:1: rule__ElkPort__PropertiesAssignment_2_2 : ( ruleProperty ) ;
    public final void rule__ElkPort__PropertiesAssignment_2_2() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:5306:1: ( ( ruleProperty ) )
            // InternalElkGraph.g:5307:2: ( ruleProperty )
            {
            // InternalElkGraph.g:5307:2: ( ruleProperty )
            // InternalElkGraph.g:5308:3: ruleProperty
            {
             before(grammarAccess.getElkPortAccess().getPropertiesPropertyParserRuleCall_2_2_0()); 
            pushFollow(FOLLOW_2);
            ruleProperty();

            state._fsp--;

             after(grammarAccess.getElkPortAccess().getPropertiesPropertyParserRuleCall_2_2_0()); 

            }


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ElkPort__PropertiesAssignment_2_2"


    // $ANTLR start "rule__ElkPort__LabelsAssignment_2_3"
    // InternalElkGraph.g:5317:1: rule__ElkPort__LabelsAssignment_2_3 : ( ruleElkLabel ) ;
    public final void rule__ElkPort__LabelsAssignment_2_3() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:5321:1: ( ( ruleElkLabel ) )
            // InternalElkGraph.g:5322:2: ( ruleElkLabel )
            {
            // InternalElkGraph.g:5322:2: ( ruleElkLabel )
            // InternalElkGraph.g:5323:3: ruleElkLabel
            {
             before(grammarAccess.getElkPortAccess().getLabelsElkLabelParserRuleCall_2_3_0()); 
            pushFollow(FOLLOW_2);
            ruleElkLabel();

            state._fsp--;

             after(grammarAccess.getElkPortAccess().getLabelsElkLabelParserRuleCall_2_3_0()); 

            }


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ElkPort__LabelsAssignment_2_3"


    // $ANTLR start "rule__ShapeLayout__XAssignment_2_0_2"
    // InternalElkGraph.g:5332:1: rule__ShapeLayout__XAssignment_2_0_2 : ( ruleNumber ) ;
    public final void rule__ShapeLayout__XAssignment_2_0_2() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:5336:1: ( ( ruleNumber ) )
            // InternalElkGraph.g:5337:2: ( ruleNumber )
            {
            // InternalElkGraph.g:5337:2: ( ruleNumber )
            // InternalElkGraph.g:5338:3: ruleNumber
            {
             before(grammarAccess.getShapeLayoutAccess().getXNumberParserRuleCall_2_0_2_0()); 
            pushFollow(FOLLOW_2);
            ruleNumber();

            state._fsp--;

             after(grammarAccess.getShapeLayoutAccess().getXNumberParserRuleCall_2_0_2_0()); 

            }


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ShapeLayout__XAssignment_2_0_2"


    // $ANTLR start "rule__ShapeLayout__YAssignment_2_0_4"
    // InternalElkGraph.g:5347:1: rule__ShapeLayout__YAssignment_2_0_4 : ( ruleNumber ) ;
    public final void rule__ShapeLayout__YAssignment_2_0_4() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:5351:1: ( ( ruleNumber ) )
            // InternalElkGraph.g:5352:2: ( ruleNumber )
            {
            // InternalElkGraph.g:5352:2: ( ruleNumber )
            // InternalElkGraph.g:5353:3: ruleNumber
            {
             before(grammarAccess.getShapeLayoutAccess().getYNumberParserRuleCall_2_0_4_0()); 
            pushFollow(FOLLOW_2);
            ruleNumber();

            state._fsp--;

             after(grammarAccess.getShapeLayoutAccess().getYNumberParserRuleCall_2_0_4_0()); 

            }


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ShapeLayout__YAssignment_2_0_4"


    // $ANTLR start "rule__ShapeLayout__WidthAssignment_2_1_2"
    // InternalElkGraph.g:5362:1: rule__ShapeLayout__WidthAssignment_2_1_2 : ( ruleNumber ) ;
    public final void rule__ShapeLayout__WidthAssignment_2_1_2() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:5366:1: ( ( ruleNumber ) )
            // InternalElkGraph.g:5367:2: ( ruleNumber )
            {
            // InternalElkGraph.g:5367:2: ( ruleNumber )
            // InternalElkGraph.g:5368:3: ruleNumber
            {
             before(grammarAccess.getShapeLayoutAccess().getWidthNumberParserRuleCall_2_1_2_0()); 
            pushFollow(FOLLOW_2);
            ruleNumber();

            state._fsp--;

             after(grammarAccess.getShapeLayoutAccess().getWidthNumberParserRuleCall_2_1_2_0()); 

            }


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ShapeLayout__WidthAssignment_2_1_2"


    // $ANTLR start "rule__ShapeLayout__HeightAssignment_2_1_4"
    // InternalElkGraph.g:5377:1: rule__ShapeLayout__HeightAssignment_2_1_4 : ( ruleNumber ) ;
    public final void rule__ShapeLayout__HeightAssignment_2_1_4() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:5381:1: ( ( ruleNumber ) )
            // InternalElkGraph.g:5382:2: ( ruleNumber )
            {
            // InternalElkGraph.g:5382:2: ( ruleNumber )
            // InternalElkGraph.g:5383:3: ruleNumber
            {
             before(grammarAccess.getShapeLayoutAccess().getHeightNumberParserRuleCall_2_1_4_0()); 
            pushFollow(FOLLOW_2);
            ruleNumber();

            state._fsp--;

             after(grammarAccess.getShapeLayoutAccess().getHeightNumberParserRuleCall_2_1_4_0()); 

            }


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ShapeLayout__HeightAssignment_2_1_4"


    // $ANTLR start "rule__ElkEdge__IdentifierAssignment_1_0"
    // InternalElkGraph.g:5392:1: rule__ElkEdge__IdentifierAssignment_1_0 : ( RULE_ID ) ;
    public final void rule__ElkEdge__IdentifierAssignment_1_0() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:5396:1: ( ( RULE_ID ) )
            // InternalElkGraph.g:5397:2: ( RULE_ID )
            {
            // InternalElkGraph.g:5397:2: ( RULE_ID )
            // InternalElkGraph.g:5398:3: RULE_ID
            {
             before(grammarAccess.getElkEdgeAccess().getIdentifierIDTerminalRuleCall_1_0_0()); 
            match(input,RULE_ID,FOLLOW_2); 
             after(grammarAccess.getElkEdgeAccess().getIdentifierIDTerminalRuleCall_1_0_0()); 

            }


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ElkEdge__IdentifierAssignment_1_0"


    // $ANTLR start "rule__ElkEdge__SourcesAssignment_2"
    // InternalElkGraph.g:5407:1: rule__ElkEdge__SourcesAssignment_2 : ( ( ruleQualifiedId ) ) ;
    public final void rule__ElkEdge__SourcesAssignment_2() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:5411:1: ( ( ( ruleQualifiedId ) ) )
            // InternalElkGraph.g:5412:2: ( ( ruleQualifiedId ) )
            {
            // InternalElkGraph.g:5412:2: ( ( ruleQualifiedId ) )
            // InternalElkGraph.g:5413:3: ( ruleQualifiedId )
            {
             before(grammarAccess.getElkEdgeAccess().getSourcesElkConnectableShapeCrossReference_2_0()); 
            // InternalElkGraph.g:5414:3: ( ruleQualifiedId )
            // InternalElkGraph.g:5415:4: ruleQualifiedId
            {
             before(grammarAccess.getElkEdgeAccess().getSourcesElkConnectableShapeQualifiedIdParserRuleCall_2_0_1()); 
            pushFollow(FOLLOW_2);
            ruleQualifiedId();

            state._fsp--;

             after(grammarAccess.getElkEdgeAccess().getSourcesElkConnectableShapeQualifiedIdParserRuleCall_2_0_1()); 

            }

             after(grammarAccess.getElkEdgeAccess().getSourcesElkConnectableShapeCrossReference_2_0()); 

            }


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ElkEdge__SourcesAssignment_2"


    // $ANTLR start "rule__ElkEdge__SourcesAssignment_3_1"
    // InternalElkGraph.g:5426:1: rule__ElkEdge__SourcesAssignment_3_1 : ( ( ruleQualifiedId ) ) ;
    public final void rule__ElkEdge__SourcesAssignment_3_1() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:5430:1: ( ( ( ruleQualifiedId ) ) )
            // InternalElkGraph.g:5431:2: ( ( ruleQualifiedId ) )
            {
            // InternalElkGraph.g:5431:2: ( ( ruleQualifiedId ) )
            // InternalElkGraph.g:5432:3: ( ruleQualifiedId )
            {
             before(grammarAccess.getElkEdgeAccess().getSourcesElkConnectableShapeCrossReference_3_1_0()); 
            // InternalElkGraph.g:5433:3: ( ruleQualifiedId )
            // InternalElkGraph.g:5434:4: ruleQualifiedId
            {
             before(grammarAccess.getElkEdgeAccess().getSourcesElkConnectableShapeQualifiedIdParserRuleCall_3_1_0_1()); 
            pushFollow(FOLLOW_2);
            ruleQualifiedId();

            state._fsp--;

             after(grammarAccess.getElkEdgeAccess().getSourcesElkConnectableShapeQualifiedIdParserRuleCall_3_1_0_1()); 

            }

             after(grammarAccess.getElkEdgeAccess().getSourcesElkConnectableShapeCrossReference_3_1_0()); 

            }


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ElkEdge__SourcesAssignment_3_1"


    // $ANTLR start "rule__ElkEdge__TargetsAssignment_5"
    // InternalElkGraph.g:5445:1: rule__ElkEdge__TargetsAssignment_5 : ( ( ruleQualifiedId ) ) ;
    public final void rule__ElkEdge__TargetsAssignment_5() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:5449:1: ( ( ( ruleQualifiedId ) ) )
            // InternalElkGraph.g:5450:2: ( ( ruleQualifiedId ) )
            {
            // InternalElkGraph.g:5450:2: ( ( ruleQualifiedId ) )
            // InternalElkGraph.g:5451:3: ( ruleQualifiedId )
            {
             before(grammarAccess.getElkEdgeAccess().getTargetsElkConnectableShapeCrossReference_5_0()); 
            // InternalElkGraph.g:5452:3: ( ruleQualifiedId )
            // InternalElkGraph.g:5453:4: ruleQualifiedId
            {
             before(grammarAccess.getElkEdgeAccess().getTargetsElkConnectableShapeQualifiedIdParserRuleCall_5_0_1()); 
            pushFollow(FOLLOW_2);
            ruleQualifiedId();

            state._fsp--;

             after(grammarAccess.getElkEdgeAccess().getTargetsElkConnectableShapeQualifiedIdParserRuleCall_5_0_1()); 

            }

             after(grammarAccess.getElkEdgeAccess().getTargetsElkConnectableShapeCrossReference_5_0()); 

            }


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ElkEdge__TargetsAssignment_5"


    // $ANTLR start "rule__ElkEdge__TargetsAssignment_6_1"
    // InternalElkGraph.g:5464:1: rule__ElkEdge__TargetsAssignment_6_1 : ( ( ruleQualifiedId ) ) ;
    public final void rule__ElkEdge__TargetsAssignment_6_1() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:5468:1: ( ( ( ruleQualifiedId ) ) )
            // InternalElkGraph.g:5469:2: ( ( ruleQualifiedId ) )
            {
            // InternalElkGraph.g:5469:2: ( ( ruleQualifiedId ) )
            // InternalElkGraph.g:5470:3: ( ruleQualifiedId )
            {
             before(grammarAccess.getElkEdgeAccess().getTargetsElkConnectableShapeCrossReference_6_1_0()); 
            // InternalElkGraph.g:5471:3: ( ruleQualifiedId )
            // InternalElkGraph.g:5472:4: ruleQualifiedId
            {
             before(grammarAccess.getElkEdgeAccess().getTargetsElkConnectableShapeQualifiedIdParserRuleCall_6_1_0_1()); 
            pushFollow(FOLLOW_2);
            ruleQualifiedId();

            state._fsp--;

             after(grammarAccess.getElkEdgeAccess().getTargetsElkConnectableShapeQualifiedIdParserRuleCall_6_1_0_1()); 

            }

             after(grammarAccess.getElkEdgeAccess().getTargetsElkConnectableShapeCrossReference_6_1_0()); 

            }


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ElkEdge__TargetsAssignment_6_1"


    // $ANTLR start "rule__ElkEdge__PropertiesAssignment_7_2"
    // InternalElkGraph.g:5483:1: rule__ElkEdge__PropertiesAssignment_7_2 : ( ruleProperty ) ;
    public final void rule__ElkEdge__PropertiesAssignment_7_2() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:5487:1: ( ( ruleProperty ) )
            // InternalElkGraph.g:5488:2: ( ruleProperty )
            {
            // InternalElkGraph.g:5488:2: ( ruleProperty )
            // InternalElkGraph.g:5489:3: ruleProperty
            {
             before(grammarAccess.getElkEdgeAccess().getPropertiesPropertyParserRuleCall_7_2_0()); 
            pushFollow(FOLLOW_2);
            ruleProperty();

            state._fsp--;

             after(grammarAccess.getElkEdgeAccess().getPropertiesPropertyParserRuleCall_7_2_0()); 

            }


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ElkEdge__PropertiesAssignment_7_2"


    // $ANTLR start "rule__ElkEdge__LabelsAssignment_7_3"
    // InternalElkGraph.g:5498:1: rule__ElkEdge__LabelsAssignment_7_3 : ( ruleElkLabel ) ;
    public final void rule__ElkEdge__LabelsAssignment_7_3() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:5502:1: ( ( ruleElkLabel ) )
            // InternalElkGraph.g:5503:2: ( ruleElkLabel )
            {
            // InternalElkGraph.g:5503:2: ( ruleElkLabel )
            // InternalElkGraph.g:5504:3: ruleElkLabel
            {
             before(grammarAccess.getElkEdgeAccess().getLabelsElkLabelParserRuleCall_7_3_0()); 
            pushFollow(FOLLOW_2);
            ruleElkLabel();

            state._fsp--;

             after(grammarAccess.getElkEdgeAccess().getLabelsElkLabelParserRuleCall_7_3_0()); 

            }


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ElkEdge__LabelsAssignment_7_3"


    // $ANTLR start "rule__EdgeLayout__SectionsAssignment_2_0"
    // InternalElkGraph.g:5513:1: rule__EdgeLayout__SectionsAssignment_2_0 : ( ruleElkSingleEdgeSection ) ;
    public final void rule__EdgeLayout__SectionsAssignment_2_0() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:5517:1: ( ( ruleElkSingleEdgeSection ) )
            // InternalElkGraph.g:5518:2: ( ruleElkSingleEdgeSection )
            {
            // InternalElkGraph.g:5518:2: ( ruleElkSingleEdgeSection )
            // InternalElkGraph.g:5519:3: ruleElkSingleEdgeSection
            {
             before(grammarAccess.getEdgeLayoutAccess().getSectionsElkSingleEdgeSectionParserRuleCall_2_0_0()); 
            pushFollow(FOLLOW_2);
            ruleElkSingleEdgeSection();

            state._fsp--;

             after(grammarAccess.getEdgeLayoutAccess().getSectionsElkSingleEdgeSectionParserRuleCall_2_0_0()); 

            }


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__EdgeLayout__SectionsAssignment_2_0"


    // $ANTLR start "rule__EdgeLayout__SectionsAssignment_2_1"
    // InternalElkGraph.g:5528:1: rule__EdgeLayout__SectionsAssignment_2_1 : ( ruleElkEdgeSection ) ;
    public final void rule__EdgeLayout__SectionsAssignment_2_1() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:5532:1: ( ( ruleElkEdgeSection ) )
            // InternalElkGraph.g:5533:2: ( ruleElkEdgeSection )
            {
            // InternalElkGraph.g:5533:2: ( ruleElkEdgeSection )
            // InternalElkGraph.g:5534:3: ruleElkEdgeSection
            {
             before(grammarAccess.getEdgeLayoutAccess().getSectionsElkEdgeSectionParserRuleCall_2_1_0()); 
            pushFollow(FOLLOW_2);
            ruleElkEdgeSection();

            state._fsp--;

             after(grammarAccess.getEdgeLayoutAccess().getSectionsElkEdgeSectionParserRuleCall_2_1_0()); 

            }


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__EdgeLayout__SectionsAssignment_2_1"


    // $ANTLR start "rule__ElkSingleEdgeSection__IncomingShapeAssignment_1_0_0_2"
    // InternalElkGraph.g:5543:1: rule__ElkSingleEdgeSection__IncomingShapeAssignment_1_0_0_2 : ( ( ruleQualifiedId ) ) ;
    public final void rule__ElkSingleEdgeSection__IncomingShapeAssignment_1_0_0_2() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:5547:1: ( ( ( ruleQualifiedId ) ) )
            // InternalElkGraph.g:5548:2: ( ( ruleQualifiedId ) )
            {
            // InternalElkGraph.g:5548:2: ( ( ruleQualifiedId ) )
            // InternalElkGraph.g:5549:3: ( ruleQualifiedId )
            {
             before(grammarAccess.getElkSingleEdgeSectionAccess().getIncomingShapeElkConnectableShapeCrossReference_1_0_0_2_0()); 
            // InternalElkGraph.g:5550:3: ( ruleQualifiedId )
            // InternalElkGraph.g:5551:4: ruleQualifiedId
            {
             before(grammarAccess.getElkSingleEdgeSectionAccess().getIncomingShapeElkConnectableShapeQualifiedIdParserRuleCall_1_0_0_2_0_1()); 
            pushFollow(FOLLOW_2);
            ruleQualifiedId();

            state._fsp--;

             after(grammarAccess.getElkSingleEdgeSectionAccess().getIncomingShapeElkConnectableShapeQualifiedIdParserRuleCall_1_0_0_2_0_1()); 

            }

             after(grammarAccess.getElkSingleEdgeSectionAccess().getIncomingShapeElkConnectableShapeCrossReference_1_0_0_2_0()); 

            }


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ElkSingleEdgeSection__IncomingShapeAssignment_1_0_0_2"


    // $ANTLR start "rule__ElkSingleEdgeSection__OutgoingShapeAssignment_1_0_1_2"
    // InternalElkGraph.g:5562:1: rule__ElkSingleEdgeSection__OutgoingShapeAssignment_1_0_1_2 : ( ( ruleQualifiedId ) ) ;
    public final void rule__ElkSingleEdgeSection__OutgoingShapeAssignment_1_0_1_2() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:5566:1: ( ( ( ruleQualifiedId ) ) )
            // InternalElkGraph.g:5567:2: ( ( ruleQualifiedId ) )
            {
            // InternalElkGraph.g:5567:2: ( ( ruleQualifiedId ) )
            // InternalElkGraph.g:5568:3: ( ruleQualifiedId )
            {
             before(grammarAccess.getElkSingleEdgeSectionAccess().getOutgoingShapeElkConnectableShapeCrossReference_1_0_1_2_0()); 
            // InternalElkGraph.g:5569:3: ( ruleQualifiedId )
            // InternalElkGraph.g:5570:4: ruleQualifiedId
            {
             before(grammarAccess.getElkSingleEdgeSectionAccess().getOutgoingShapeElkConnectableShapeQualifiedIdParserRuleCall_1_0_1_2_0_1()); 
            pushFollow(FOLLOW_2);
            ruleQualifiedId();

            state._fsp--;

             after(grammarAccess.getElkSingleEdgeSectionAccess().getOutgoingShapeElkConnectableShapeQualifiedIdParserRuleCall_1_0_1_2_0_1()); 

            }

             after(grammarAccess.getElkSingleEdgeSectionAccess().getOutgoingShapeElkConnectableShapeCrossReference_1_0_1_2_0()); 

            }


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ElkSingleEdgeSection__OutgoingShapeAssignment_1_0_1_2"


    // $ANTLR start "rule__ElkSingleEdgeSection__StartXAssignment_1_0_2_2"
    // InternalElkGraph.g:5581:1: rule__ElkSingleEdgeSection__StartXAssignment_1_0_2_2 : ( ruleNumber ) ;
    public final void rule__ElkSingleEdgeSection__StartXAssignment_1_0_2_2() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:5585:1: ( ( ruleNumber ) )
            // InternalElkGraph.g:5586:2: ( ruleNumber )
            {
            // InternalElkGraph.g:5586:2: ( ruleNumber )
            // InternalElkGraph.g:5587:3: ruleNumber
            {
             before(grammarAccess.getElkSingleEdgeSectionAccess().getStartXNumberParserRuleCall_1_0_2_2_0()); 
            pushFollow(FOLLOW_2);
            ruleNumber();

            state._fsp--;

             after(grammarAccess.getElkSingleEdgeSectionAccess().getStartXNumberParserRuleCall_1_0_2_2_0()); 

            }


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ElkSingleEdgeSection__StartXAssignment_1_0_2_2"


    // $ANTLR start "rule__ElkSingleEdgeSection__StartYAssignment_1_0_2_4"
    // InternalElkGraph.g:5596:1: rule__ElkSingleEdgeSection__StartYAssignment_1_0_2_4 : ( ruleNumber ) ;
    public final void rule__ElkSingleEdgeSection__StartYAssignment_1_0_2_4() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:5600:1: ( ( ruleNumber ) )
            // InternalElkGraph.g:5601:2: ( ruleNumber )
            {
            // InternalElkGraph.g:5601:2: ( ruleNumber )
            // InternalElkGraph.g:5602:3: ruleNumber
            {
             before(grammarAccess.getElkSingleEdgeSectionAccess().getStartYNumberParserRuleCall_1_0_2_4_0()); 
            pushFollow(FOLLOW_2);
            ruleNumber();

            state._fsp--;

             after(grammarAccess.getElkSingleEdgeSectionAccess().getStartYNumberParserRuleCall_1_0_2_4_0()); 

            }


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ElkSingleEdgeSection__StartYAssignment_1_0_2_4"


    // $ANTLR start "rule__ElkSingleEdgeSection__EndXAssignment_1_0_3_2"
    // InternalElkGraph.g:5611:1: rule__ElkSingleEdgeSection__EndXAssignment_1_0_3_2 : ( ruleNumber ) ;
    public final void rule__ElkSingleEdgeSection__EndXAssignment_1_0_3_2() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:5615:1: ( ( ruleNumber ) )
            // InternalElkGraph.g:5616:2: ( ruleNumber )
            {
            // InternalElkGraph.g:5616:2: ( ruleNumber )
            // InternalElkGraph.g:5617:3: ruleNumber
            {
             before(grammarAccess.getElkSingleEdgeSectionAccess().getEndXNumberParserRuleCall_1_0_3_2_0()); 
            pushFollow(FOLLOW_2);
            ruleNumber();

            state._fsp--;

             after(grammarAccess.getElkSingleEdgeSectionAccess().getEndXNumberParserRuleCall_1_0_3_2_0()); 

            }


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ElkSingleEdgeSection__EndXAssignment_1_0_3_2"


    // $ANTLR start "rule__ElkSingleEdgeSection__EndYAssignment_1_0_3_4"
    // InternalElkGraph.g:5626:1: rule__ElkSingleEdgeSection__EndYAssignment_1_0_3_4 : ( ruleNumber ) ;
    public final void rule__ElkSingleEdgeSection__EndYAssignment_1_0_3_4() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:5630:1: ( ( ruleNumber ) )
            // InternalElkGraph.g:5631:2: ( ruleNumber )
            {
            // InternalElkGraph.g:5631:2: ( ruleNumber )
            // InternalElkGraph.g:5632:3: ruleNumber
            {
             before(grammarAccess.getElkSingleEdgeSectionAccess().getEndYNumberParserRuleCall_1_0_3_4_0()); 
            pushFollow(FOLLOW_2);
            ruleNumber();

            state._fsp--;

             after(grammarAccess.getElkSingleEdgeSectionAccess().getEndYNumberParserRuleCall_1_0_3_4_0()); 

            }


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ElkSingleEdgeSection__EndYAssignment_1_0_3_4"


    // $ANTLR start "rule__ElkSingleEdgeSection__BendPointsAssignment_1_1_2"
    // InternalElkGraph.g:5641:1: rule__ElkSingleEdgeSection__BendPointsAssignment_1_1_2 : ( ruleElkBendPoint ) ;
    public final void rule__ElkSingleEdgeSection__BendPointsAssignment_1_1_2() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:5645:1: ( ( ruleElkBendPoint ) )
            // InternalElkGraph.g:5646:2: ( ruleElkBendPoint )
            {
            // InternalElkGraph.g:5646:2: ( ruleElkBendPoint )
            // InternalElkGraph.g:5647:3: ruleElkBendPoint
            {
             before(grammarAccess.getElkSingleEdgeSectionAccess().getBendPointsElkBendPointParserRuleCall_1_1_2_0()); 
            pushFollow(FOLLOW_2);
            ruleElkBendPoint();

            state._fsp--;

             after(grammarAccess.getElkSingleEdgeSectionAccess().getBendPointsElkBendPointParserRuleCall_1_1_2_0()); 

            }


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ElkSingleEdgeSection__BendPointsAssignment_1_1_2"


    // $ANTLR start "rule__ElkSingleEdgeSection__BendPointsAssignment_1_1_3_1"
    // InternalElkGraph.g:5656:1: rule__ElkSingleEdgeSection__BendPointsAssignment_1_1_3_1 : ( ruleElkBendPoint ) ;
    public final void rule__ElkSingleEdgeSection__BendPointsAssignment_1_1_3_1() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:5660:1: ( ( ruleElkBendPoint ) )
            // InternalElkGraph.g:5661:2: ( ruleElkBendPoint )
            {
            // InternalElkGraph.g:5661:2: ( ruleElkBendPoint )
            // InternalElkGraph.g:5662:3: ruleElkBendPoint
            {
             before(grammarAccess.getElkSingleEdgeSectionAccess().getBendPointsElkBendPointParserRuleCall_1_1_3_1_0()); 
            pushFollow(FOLLOW_2);
            ruleElkBendPoint();

            state._fsp--;

             after(grammarAccess.getElkSingleEdgeSectionAccess().getBendPointsElkBendPointParserRuleCall_1_1_3_1_0()); 

            }


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ElkSingleEdgeSection__BendPointsAssignment_1_1_3_1"


    // $ANTLR start "rule__ElkSingleEdgeSection__PropertiesAssignment_1_2"
    // InternalElkGraph.g:5671:1: rule__ElkSingleEdgeSection__PropertiesAssignment_1_2 : ( ruleProperty ) ;
    public final void rule__ElkSingleEdgeSection__PropertiesAssignment_1_2() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:5675:1: ( ( ruleProperty ) )
            // InternalElkGraph.g:5676:2: ( ruleProperty )
            {
            // InternalElkGraph.g:5676:2: ( ruleProperty )
            // InternalElkGraph.g:5677:3: ruleProperty
            {
             before(grammarAccess.getElkSingleEdgeSectionAccess().getPropertiesPropertyParserRuleCall_1_2_0()); 
            pushFollow(FOLLOW_2);
            ruleProperty();

            state._fsp--;

             after(grammarAccess.getElkSingleEdgeSectionAccess().getPropertiesPropertyParserRuleCall_1_2_0()); 

            }


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ElkSingleEdgeSection__PropertiesAssignment_1_2"


    // $ANTLR start "rule__ElkEdgeSection__IdentifierAssignment_1"
    // InternalElkGraph.g:5686:1: rule__ElkEdgeSection__IdentifierAssignment_1 : ( RULE_ID ) ;
    public final void rule__ElkEdgeSection__IdentifierAssignment_1() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:5690:1: ( ( RULE_ID ) )
            // InternalElkGraph.g:5691:2: ( RULE_ID )
            {
            // InternalElkGraph.g:5691:2: ( RULE_ID )
            // InternalElkGraph.g:5692:3: RULE_ID
            {
             before(grammarAccess.getElkEdgeSectionAccess().getIdentifierIDTerminalRuleCall_1_0()); 
            match(input,RULE_ID,FOLLOW_2); 
             after(grammarAccess.getElkEdgeSectionAccess().getIdentifierIDTerminalRuleCall_1_0()); 

            }


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ElkEdgeSection__IdentifierAssignment_1"


    // $ANTLR start "rule__ElkEdgeSection__OutgoingSectionsAssignment_2_1"
    // InternalElkGraph.g:5701:1: rule__ElkEdgeSection__OutgoingSectionsAssignment_2_1 : ( ( RULE_ID ) ) ;
    public final void rule__ElkEdgeSection__OutgoingSectionsAssignment_2_1() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:5705:1: ( ( ( RULE_ID ) ) )
            // InternalElkGraph.g:5706:2: ( ( RULE_ID ) )
            {
            // InternalElkGraph.g:5706:2: ( ( RULE_ID ) )
            // InternalElkGraph.g:5707:3: ( RULE_ID )
            {
             before(grammarAccess.getElkEdgeSectionAccess().getOutgoingSectionsElkEdgeSectionCrossReference_2_1_0()); 
            // InternalElkGraph.g:5708:3: ( RULE_ID )
            // InternalElkGraph.g:5709:4: RULE_ID
            {
             before(grammarAccess.getElkEdgeSectionAccess().getOutgoingSectionsElkEdgeSectionIDTerminalRuleCall_2_1_0_1()); 
            match(input,RULE_ID,FOLLOW_2); 
             after(grammarAccess.getElkEdgeSectionAccess().getOutgoingSectionsElkEdgeSectionIDTerminalRuleCall_2_1_0_1()); 

            }

             after(grammarAccess.getElkEdgeSectionAccess().getOutgoingSectionsElkEdgeSectionCrossReference_2_1_0()); 

            }


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ElkEdgeSection__OutgoingSectionsAssignment_2_1"


    // $ANTLR start "rule__ElkEdgeSection__OutgoingSectionsAssignment_2_2_1"
    // InternalElkGraph.g:5720:1: rule__ElkEdgeSection__OutgoingSectionsAssignment_2_2_1 : ( ( RULE_ID ) ) ;
    public final void rule__ElkEdgeSection__OutgoingSectionsAssignment_2_2_1() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:5724:1: ( ( ( RULE_ID ) ) )
            // InternalElkGraph.g:5725:2: ( ( RULE_ID ) )
            {
            // InternalElkGraph.g:5725:2: ( ( RULE_ID ) )
            // InternalElkGraph.g:5726:3: ( RULE_ID )
            {
             before(grammarAccess.getElkEdgeSectionAccess().getOutgoingSectionsElkEdgeSectionCrossReference_2_2_1_0()); 
            // InternalElkGraph.g:5727:3: ( RULE_ID )
            // InternalElkGraph.g:5728:4: RULE_ID
            {
             before(grammarAccess.getElkEdgeSectionAccess().getOutgoingSectionsElkEdgeSectionIDTerminalRuleCall_2_2_1_0_1()); 
            match(input,RULE_ID,FOLLOW_2); 
             after(grammarAccess.getElkEdgeSectionAccess().getOutgoingSectionsElkEdgeSectionIDTerminalRuleCall_2_2_1_0_1()); 

            }

             after(grammarAccess.getElkEdgeSectionAccess().getOutgoingSectionsElkEdgeSectionCrossReference_2_2_1_0()); 

            }


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ElkEdgeSection__OutgoingSectionsAssignment_2_2_1"


    // $ANTLR start "rule__ElkEdgeSection__IncomingShapeAssignment_4_0_0_2"
    // InternalElkGraph.g:5739:1: rule__ElkEdgeSection__IncomingShapeAssignment_4_0_0_2 : ( ( ruleQualifiedId ) ) ;
    public final void rule__ElkEdgeSection__IncomingShapeAssignment_4_0_0_2() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:5743:1: ( ( ( ruleQualifiedId ) ) )
            // InternalElkGraph.g:5744:2: ( ( ruleQualifiedId ) )
            {
            // InternalElkGraph.g:5744:2: ( ( ruleQualifiedId ) )
            // InternalElkGraph.g:5745:3: ( ruleQualifiedId )
            {
             before(grammarAccess.getElkEdgeSectionAccess().getIncomingShapeElkConnectableShapeCrossReference_4_0_0_2_0()); 
            // InternalElkGraph.g:5746:3: ( ruleQualifiedId )
            // InternalElkGraph.g:5747:4: ruleQualifiedId
            {
             before(grammarAccess.getElkEdgeSectionAccess().getIncomingShapeElkConnectableShapeQualifiedIdParserRuleCall_4_0_0_2_0_1()); 
            pushFollow(FOLLOW_2);
            ruleQualifiedId();

            state._fsp--;

             after(grammarAccess.getElkEdgeSectionAccess().getIncomingShapeElkConnectableShapeQualifiedIdParserRuleCall_4_0_0_2_0_1()); 

            }

             after(grammarAccess.getElkEdgeSectionAccess().getIncomingShapeElkConnectableShapeCrossReference_4_0_0_2_0()); 

            }


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ElkEdgeSection__IncomingShapeAssignment_4_0_0_2"


    // $ANTLR start "rule__ElkEdgeSection__OutgoingShapeAssignment_4_0_1_2"
    // InternalElkGraph.g:5758:1: rule__ElkEdgeSection__OutgoingShapeAssignment_4_0_1_2 : ( ( ruleQualifiedId ) ) ;
    public final void rule__ElkEdgeSection__OutgoingShapeAssignment_4_0_1_2() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:5762:1: ( ( ( ruleQualifiedId ) ) )
            // InternalElkGraph.g:5763:2: ( ( ruleQualifiedId ) )
            {
            // InternalElkGraph.g:5763:2: ( ( ruleQualifiedId ) )
            // InternalElkGraph.g:5764:3: ( ruleQualifiedId )
            {
             before(grammarAccess.getElkEdgeSectionAccess().getOutgoingShapeElkConnectableShapeCrossReference_4_0_1_2_0()); 
            // InternalElkGraph.g:5765:3: ( ruleQualifiedId )
            // InternalElkGraph.g:5766:4: ruleQualifiedId
            {
             before(grammarAccess.getElkEdgeSectionAccess().getOutgoingShapeElkConnectableShapeQualifiedIdParserRuleCall_4_0_1_2_0_1()); 
            pushFollow(FOLLOW_2);
            ruleQualifiedId();

            state._fsp--;

             after(grammarAccess.getElkEdgeSectionAccess().getOutgoingShapeElkConnectableShapeQualifiedIdParserRuleCall_4_0_1_2_0_1()); 

            }

             after(grammarAccess.getElkEdgeSectionAccess().getOutgoingShapeElkConnectableShapeCrossReference_4_0_1_2_0()); 

            }


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ElkEdgeSection__OutgoingShapeAssignment_4_0_1_2"


    // $ANTLR start "rule__ElkEdgeSection__StartXAssignment_4_0_2_2"
    // InternalElkGraph.g:5777:1: rule__ElkEdgeSection__StartXAssignment_4_0_2_2 : ( ruleNumber ) ;
    public final void rule__ElkEdgeSection__StartXAssignment_4_0_2_2() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:5781:1: ( ( ruleNumber ) )
            // InternalElkGraph.g:5782:2: ( ruleNumber )
            {
            // InternalElkGraph.g:5782:2: ( ruleNumber )
            // InternalElkGraph.g:5783:3: ruleNumber
            {
             before(grammarAccess.getElkEdgeSectionAccess().getStartXNumberParserRuleCall_4_0_2_2_0()); 
            pushFollow(FOLLOW_2);
            ruleNumber();

            state._fsp--;

             after(grammarAccess.getElkEdgeSectionAccess().getStartXNumberParserRuleCall_4_0_2_2_0()); 

            }


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ElkEdgeSection__StartXAssignment_4_0_2_2"


    // $ANTLR start "rule__ElkEdgeSection__StartYAssignment_4_0_2_4"
    // InternalElkGraph.g:5792:1: rule__ElkEdgeSection__StartYAssignment_4_0_2_4 : ( ruleNumber ) ;
    public final void rule__ElkEdgeSection__StartYAssignment_4_0_2_4() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:5796:1: ( ( ruleNumber ) )
            // InternalElkGraph.g:5797:2: ( ruleNumber )
            {
            // InternalElkGraph.g:5797:2: ( ruleNumber )
            // InternalElkGraph.g:5798:3: ruleNumber
            {
             before(grammarAccess.getElkEdgeSectionAccess().getStartYNumberParserRuleCall_4_0_2_4_0()); 
            pushFollow(FOLLOW_2);
            ruleNumber();

            state._fsp--;

             after(grammarAccess.getElkEdgeSectionAccess().getStartYNumberParserRuleCall_4_0_2_4_0()); 

            }


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ElkEdgeSection__StartYAssignment_4_0_2_4"


    // $ANTLR start "rule__ElkEdgeSection__EndXAssignment_4_0_3_2"
    // InternalElkGraph.g:5807:1: rule__ElkEdgeSection__EndXAssignment_4_0_3_2 : ( ruleNumber ) ;
    public final void rule__ElkEdgeSection__EndXAssignment_4_0_3_2() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:5811:1: ( ( ruleNumber ) )
            // InternalElkGraph.g:5812:2: ( ruleNumber )
            {
            // InternalElkGraph.g:5812:2: ( ruleNumber )
            // InternalElkGraph.g:5813:3: ruleNumber
            {
             before(grammarAccess.getElkEdgeSectionAccess().getEndXNumberParserRuleCall_4_0_3_2_0()); 
            pushFollow(FOLLOW_2);
            ruleNumber();

            state._fsp--;

             after(grammarAccess.getElkEdgeSectionAccess().getEndXNumberParserRuleCall_4_0_3_2_0()); 

            }


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ElkEdgeSection__EndXAssignment_4_0_3_2"


    // $ANTLR start "rule__ElkEdgeSection__EndYAssignment_4_0_3_4"
    // InternalElkGraph.g:5822:1: rule__ElkEdgeSection__EndYAssignment_4_0_3_4 : ( ruleNumber ) ;
    public final void rule__ElkEdgeSection__EndYAssignment_4_0_3_4() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:5826:1: ( ( ruleNumber ) )
            // InternalElkGraph.g:5827:2: ( ruleNumber )
            {
            // InternalElkGraph.g:5827:2: ( ruleNumber )
            // InternalElkGraph.g:5828:3: ruleNumber
            {
             before(grammarAccess.getElkEdgeSectionAccess().getEndYNumberParserRuleCall_4_0_3_4_0()); 
            pushFollow(FOLLOW_2);
            ruleNumber();

            state._fsp--;

             after(grammarAccess.getElkEdgeSectionAccess().getEndYNumberParserRuleCall_4_0_3_4_0()); 

            }


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ElkEdgeSection__EndYAssignment_4_0_3_4"


    // $ANTLR start "rule__ElkEdgeSection__BendPointsAssignment_4_1_2"
    // InternalElkGraph.g:5837:1: rule__ElkEdgeSection__BendPointsAssignment_4_1_2 : ( ruleElkBendPoint ) ;
    public final void rule__ElkEdgeSection__BendPointsAssignment_4_1_2() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:5841:1: ( ( ruleElkBendPoint ) )
            // InternalElkGraph.g:5842:2: ( ruleElkBendPoint )
            {
            // InternalElkGraph.g:5842:2: ( ruleElkBendPoint )
            // InternalElkGraph.g:5843:3: ruleElkBendPoint
            {
             before(grammarAccess.getElkEdgeSectionAccess().getBendPointsElkBendPointParserRuleCall_4_1_2_0()); 
            pushFollow(FOLLOW_2);
            ruleElkBendPoint();

            state._fsp--;

             after(grammarAccess.getElkEdgeSectionAccess().getBendPointsElkBendPointParserRuleCall_4_1_2_0()); 

            }


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ElkEdgeSection__BendPointsAssignment_4_1_2"


    // $ANTLR start "rule__ElkEdgeSection__BendPointsAssignment_4_1_3_1"
    // InternalElkGraph.g:5852:1: rule__ElkEdgeSection__BendPointsAssignment_4_1_3_1 : ( ruleElkBendPoint ) ;
    public final void rule__ElkEdgeSection__BendPointsAssignment_4_1_3_1() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:5856:1: ( ( ruleElkBendPoint ) )
            // InternalElkGraph.g:5857:2: ( ruleElkBendPoint )
            {
            // InternalElkGraph.g:5857:2: ( ruleElkBendPoint )
            // InternalElkGraph.g:5858:3: ruleElkBendPoint
            {
             before(grammarAccess.getElkEdgeSectionAccess().getBendPointsElkBendPointParserRuleCall_4_1_3_1_0()); 
            pushFollow(FOLLOW_2);
            ruleElkBendPoint();

            state._fsp--;

             after(grammarAccess.getElkEdgeSectionAccess().getBendPointsElkBendPointParserRuleCall_4_1_3_1_0()); 

            }


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ElkEdgeSection__BendPointsAssignment_4_1_3_1"


    // $ANTLR start "rule__ElkEdgeSection__PropertiesAssignment_4_2"
    // InternalElkGraph.g:5867:1: rule__ElkEdgeSection__PropertiesAssignment_4_2 : ( ruleProperty ) ;
    public final void rule__ElkEdgeSection__PropertiesAssignment_4_2() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:5871:1: ( ( ruleProperty ) )
            // InternalElkGraph.g:5872:2: ( ruleProperty )
            {
            // InternalElkGraph.g:5872:2: ( ruleProperty )
            // InternalElkGraph.g:5873:3: ruleProperty
            {
             before(grammarAccess.getElkEdgeSectionAccess().getPropertiesPropertyParserRuleCall_4_2_0()); 
            pushFollow(FOLLOW_2);
            ruleProperty();

            state._fsp--;

             after(grammarAccess.getElkEdgeSectionAccess().getPropertiesPropertyParserRuleCall_4_2_0()); 

            }


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ElkEdgeSection__PropertiesAssignment_4_2"


    // $ANTLR start "rule__ElkBendPoint__XAssignment_0"
    // InternalElkGraph.g:5882:1: rule__ElkBendPoint__XAssignment_0 : ( ruleNumber ) ;
    public final void rule__ElkBendPoint__XAssignment_0() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:5886:1: ( ( ruleNumber ) )
            // InternalElkGraph.g:5887:2: ( ruleNumber )
            {
            // InternalElkGraph.g:5887:2: ( ruleNumber )
            // InternalElkGraph.g:5888:3: ruleNumber
            {
             before(grammarAccess.getElkBendPointAccess().getXNumberParserRuleCall_0_0()); 
            pushFollow(FOLLOW_2);
            ruleNumber();

            state._fsp--;

             after(grammarAccess.getElkBendPointAccess().getXNumberParserRuleCall_0_0()); 

            }


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ElkBendPoint__XAssignment_0"


    // $ANTLR start "rule__ElkBendPoint__YAssignment_2"
    // InternalElkGraph.g:5897:1: rule__ElkBendPoint__YAssignment_2 : ( ruleNumber ) ;
    public final void rule__ElkBendPoint__YAssignment_2() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:5901:1: ( ( ruleNumber ) )
            // InternalElkGraph.g:5902:2: ( ruleNumber )
            {
            // InternalElkGraph.g:5902:2: ( ruleNumber )
            // InternalElkGraph.g:5903:3: ruleNumber
            {
             before(grammarAccess.getElkBendPointAccess().getYNumberParserRuleCall_2_0()); 
            pushFollow(FOLLOW_2);
            ruleNumber();

            state._fsp--;

             after(grammarAccess.getElkBendPointAccess().getYNumberParserRuleCall_2_0()); 

            }


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ElkBendPoint__YAssignment_2"


    // $ANTLR start "rule__Property__KeyAssignment_0"
    // InternalElkGraph.g:5912:1: rule__Property__KeyAssignment_0 : ( rulePropertyKey ) ;
    public final void rule__Property__KeyAssignment_0() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:5916:1: ( ( rulePropertyKey ) )
            // InternalElkGraph.g:5917:2: ( rulePropertyKey )
            {
            // InternalElkGraph.g:5917:2: ( rulePropertyKey )
            // InternalElkGraph.g:5918:3: rulePropertyKey
            {
             before(grammarAccess.getPropertyAccess().getKeyPropertyKeyParserRuleCall_0_0()); 
            pushFollow(FOLLOW_2);
            rulePropertyKey();

            state._fsp--;

             after(grammarAccess.getPropertyAccess().getKeyPropertyKeyParserRuleCall_0_0()); 

            }


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__Property__KeyAssignment_0"


    // $ANTLR start "rule__Property__ValueAssignment_2_0"
    // InternalElkGraph.g:5927:1: rule__Property__ValueAssignment_2_0 : ( ruleStringValue ) ;
    public final void rule__Property__ValueAssignment_2_0() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:5931:1: ( ( ruleStringValue ) )
            // InternalElkGraph.g:5932:2: ( ruleStringValue )
            {
            // InternalElkGraph.g:5932:2: ( ruleStringValue )
            // InternalElkGraph.g:5933:3: ruleStringValue
            {
             before(grammarAccess.getPropertyAccess().getValueStringValueParserRuleCall_2_0_0()); 
            pushFollow(FOLLOW_2);
            ruleStringValue();

            state._fsp--;

             after(grammarAccess.getPropertyAccess().getValueStringValueParserRuleCall_2_0_0()); 

            }


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__Property__ValueAssignment_2_0"


    // $ANTLR start "rule__Property__ValueAssignment_2_1"
    // InternalElkGraph.g:5942:1: rule__Property__ValueAssignment_2_1 : ( ruleQualifiedIdValue ) ;
    public final void rule__Property__ValueAssignment_2_1() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:5946:1: ( ( ruleQualifiedIdValue ) )
            // InternalElkGraph.g:5947:2: ( ruleQualifiedIdValue )
            {
            // InternalElkGraph.g:5947:2: ( ruleQualifiedIdValue )
            // InternalElkGraph.g:5948:3: ruleQualifiedIdValue
            {
             before(grammarAccess.getPropertyAccess().getValueQualifiedIdValueParserRuleCall_2_1_0()); 
            pushFollow(FOLLOW_2);
            ruleQualifiedIdValue();

            state._fsp--;

             after(grammarAccess.getPropertyAccess().getValueQualifiedIdValueParserRuleCall_2_1_0()); 

            }


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__Property__ValueAssignment_2_1"


    // $ANTLR start "rule__Property__ValueAssignment_2_2"
    // InternalElkGraph.g:5957:1: rule__Property__ValueAssignment_2_2 : ( ruleNumberValue ) ;
    public final void rule__Property__ValueAssignment_2_2() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:5961:1: ( ( ruleNumberValue ) )
            // InternalElkGraph.g:5962:2: ( ruleNumberValue )
            {
            // InternalElkGraph.g:5962:2: ( ruleNumberValue )
            // InternalElkGraph.g:5963:3: ruleNumberValue
            {
             before(grammarAccess.getPropertyAccess().getValueNumberValueParserRuleCall_2_2_0()); 
            pushFollow(FOLLOW_2);
            ruleNumberValue();

            state._fsp--;

             after(grammarAccess.getPropertyAccess().getValueNumberValueParserRuleCall_2_2_0()); 

            }


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__Property__ValueAssignment_2_2"


    // $ANTLR start "rule__Property__ValueAssignment_2_3"
    // InternalElkGraph.g:5972:1: rule__Property__ValueAssignment_2_3 : ( ruleBooleanValue ) ;
    public final void rule__Property__ValueAssignment_2_3() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraph.g:5976:1: ( ( ruleBooleanValue ) )
            // InternalElkGraph.g:5977:2: ( ruleBooleanValue )
            {
            // InternalElkGraph.g:5977:2: ( ruleBooleanValue )
            // InternalElkGraph.g:5978:3: ruleBooleanValue
            {
             before(grammarAccess.getPropertyAccess().getValueBooleanValueParserRuleCall_2_3_0()); 
            pushFollow(FOLLOW_2);
            ruleBooleanValue();

            state._fsp--;

             after(grammarAccess.getPropertyAccess().getValueBooleanValueParserRuleCall_2_3_0()); 

            }


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__Property__ValueAssignment_2_3"

    // Delegated rules


 

    public static final BitSet FOLLOW_1 = new BitSet(new long[]{0x0000000000000000L});
    public static final BitSet FOLLOW_2 = new BitSet(new long[]{0x0000000000000002L});
    public static final BitSet FOLLOW_3 = new BitSet(new long[]{0x0000002000000002L});
    public static final BitSet FOLLOW_4 = new BitSet(new long[]{0x0000000020D30080L});
    public static final BitSet FOLLOW_5 = new BitSet(new long[]{0x0000000000000082L});
    public static final BitSet FOLLOW_6 = new BitSet(new long[]{0x0000000020520002L});
    public static final BitSet FOLLOW_7 = new BitSet(new long[]{0x0000000000000080L});
    public static final BitSet FOLLOW_8 = new BitSet(new long[]{0x0000000000040000L});
    public static final BitSet FOLLOW_9 = new BitSet(new long[]{0x0000000020DA0080L});
    public static final BitSet FOLLOW_10 = new BitSet(new long[]{0x0000000000000090L});
    public static final BitSet FOLLOW_11 = new BitSet(new long[]{0x0000000000200000L});
    public static final BitSet FOLLOW_12 = new BitSet(new long[]{0x0000000000980080L});
    public static final BitSet FOLLOW_13 = new BitSet(new long[]{0x0000000000100002L});
    public static final BitSet FOLLOW_14 = new BitSet(new long[]{0x0000000001000000L});
    public static final BitSet FOLLOW_15 = new BitSet(new long[]{0x0000000014000000L});
    public static final BitSet FOLLOW_16 = new BitSet(new long[]{0x0000000002000000L});
    public static final BitSet FOLLOW_17 = new BitSet(new long[]{0x0000000000000060L});
    public static final BitSet FOLLOW_18 = new BitSet(new long[]{0x0000000008000000L});
    public static final BitSet FOLLOW_19 = new BitSet(new long[]{0x0000000048000000L});
    public static final BitSet FOLLOW_20 = new BitSet(new long[]{0x0000000008000002L});
    public static final BitSet FOLLOW_21 = new BitSet(new long[]{0x0000000008040000L});
    public static final BitSet FOLLOW_22 = new BitSet(new long[]{0x0000002780000000L});
    public static final BitSet FOLLOW_23 = new BitSet(new long[]{0x0000000780000000L});
    public static final BitSet FOLLOW_24 = new BitSet(new long[]{0x0000000800000080L});
    public static final BitSet FOLLOW_25 = new BitSet(new long[]{0x0000001000000000L});
    public static final BitSet FOLLOW_26 = new BitSet(new long[]{0x0000001000000002L});
    public static final BitSet FOLLOW_27 = new BitSet(new long[]{0x0000000041000000L});
    public static final BitSet FOLLOW_28 = new BitSet(new long[]{0x0000004000000000L});
    public static final BitSet FOLLOW_29 = new BitSet(new long[]{0x0000004000000002L});
    public static final BitSet FOLLOW_30 = new BitSet(new long[]{0x000000000000E0F0L});
    public static final BitSet FOLLOW_31 = new BitSet(new long[]{0x0000000014000002L});
    public static final BitSet FOLLOW_32 = new BitSet(new long[]{0x0000000780000002L});

}
