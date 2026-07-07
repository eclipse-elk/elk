package org.eclipse.elk.graph.json.text.ide.contentassist.antlr.internal;

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
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
@SuppressWarnings("all")
public class InternalElkGraphJsonParser extends AbstractInternalContentAssistParser {
    public static final String[] tokenNames = new String[] {
        "<invalid>", "<EOR>", "<DOWN>", "<UP>", "RULE_STRING", "RULE_ID", "RULE_SIGNED_INT", "RULE_FLOAT", "RULE_INT", "RULE_ML_COMMENT", "RULE_SL_COMMENT", "RULE_WS", "RULE_ANY_OTHER", "'true'", "'false'", "'null'", "'\"children\"'", "'\\'children\\''", "'children'", "'\"ports\"'", "'\\'ports\\''", "'ports'", "'\"labels\"'", "'\\'labels\\''", "'labels'", "'\"edges\"'", "'\\'edges\\''", "'edges'", "'\"layoutOptions\"'", "'\\'layoutOptions\\''", "'layoutOptions'", "'\"properties\"'", "'\\'properties\\''", "'properties'", "'\"id\"'", "'\\'id\\''", "'id'", "'\"x\"'", "'\\'x\\''", "'x'", "'\"y\"'", "'\\'y\\''", "'y'", "'\"width\"'", "'\\'width\\''", "'width'", "'\"height\"'", "'\\'height\\''", "'height'", "'\"sources\"'", "'\\'sources\\''", "'sources'", "'\"targets\"'", "'\\'targets\\''", "'targets'", "'\"text\"'", "'\\'text\\''", "'text'", "'{'", "','", "'}'", "':'", "'['", "']'"
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

    	public void setGrammarAccess(ElkGraphJsonGrammarAccess grammarAccess) {
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



    // $ANTLR start "entryRuleElkNode"
    // InternalElkGraphJson.g:59:1: entryRuleElkNode : ruleElkNode EOF ;
    public final void entryRuleElkNode() throws RecognitionException {
        try {
            // InternalElkGraphJson.g:60:1: ( ruleElkNode EOF )
            // InternalElkGraphJson.g:61:1: ruleElkNode EOF
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
    // InternalElkGraphJson.g:68:1: ruleElkNode : ( ( rule__ElkNode__Group__0 ) ) ;
    public final void ruleElkNode() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:72:2: ( ( ( rule__ElkNode__Group__0 ) ) )
            // InternalElkGraphJson.g:73:2: ( ( rule__ElkNode__Group__0 ) )
            {
            // InternalElkGraphJson.g:73:2: ( ( rule__ElkNode__Group__0 ) )
            // InternalElkGraphJson.g:74:3: ( rule__ElkNode__Group__0 )
            {
             before(grammarAccess.getElkNodeAccess().getGroup()); 
            // InternalElkGraphJson.g:75:3: ( rule__ElkNode__Group__0 )
            // InternalElkGraphJson.g:75:4: rule__ElkNode__Group__0
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


    // $ANTLR start "ruleNodeElement"
    // InternalElkGraphJson.g:85:1: ruleNodeElement : ( ( rule__NodeElement__Alternatives ) ) ;
    public final void ruleNodeElement() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:89:2: ( ( ( rule__NodeElement__Alternatives ) ) )
            // InternalElkGraphJson.g:90:2: ( ( rule__NodeElement__Alternatives ) )
            {
            // InternalElkGraphJson.g:90:2: ( ( rule__NodeElement__Alternatives ) )
            // InternalElkGraphJson.g:91:3: ( rule__NodeElement__Alternatives )
            {
             before(grammarAccess.getNodeElementAccess().getAlternatives()); 
            // InternalElkGraphJson.g:92:3: ( rule__NodeElement__Alternatives )
            // InternalElkGraphJson.g:92:4: rule__NodeElement__Alternatives
            {
            pushFollow(FOLLOW_2);
            rule__NodeElement__Alternatives();

            state._fsp--;


            }

             after(grammarAccess.getNodeElementAccess().getAlternatives()); 

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
    // $ANTLR end "ruleNodeElement"


    // $ANTLR start "entryRuleElkPort"
    // InternalElkGraphJson.g:101:1: entryRuleElkPort : ruleElkPort EOF ;
    public final void entryRuleElkPort() throws RecognitionException {
        try {
            // InternalElkGraphJson.g:102:1: ( ruleElkPort EOF )
            // InternalElkGraphJson.g:103:1: ruleElkPort EOF
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
    // InternalElkGraphJson.g:110:1: ruleElkPort : ( ( rule__ElkPort__Group__0 ) ) ;
    public final void ruleElkPort() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:114:2: ( ( ( rule__ElkPort__Group__0 ) ) )
            // InternalElkGraphJson.g:115:2: ( ( rule__ElkPort__Group__0 ) )
            {
            // InternalElkGraphJson.g:115:2: ( ( rule__ElkPort__Group__0 ) )
            // InternalElkGraphJson.g:116:3: ( rule__ElkPort__Group__0 )
            {
             before(grammarAccess.getElkPortAccess().getGroup()); 
            // InternalElkGraphJson.g:117:3: ( rule__ElkPort__Group__0 )
            // InternalElkGraphJson.g:117:4: rule__ElkPort__Group__0
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


    // $ANTLR start "rulePortElement"
    // InternalElkGraphJson.g:127:1: rulePortElement : ( ( rule__PortElement__Alternatives ) ) ;
    public final void rulePortElement() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:131:2: ( ( ( rule__PortElement__Alternatives ) ) )
            // InternalElkGraphJson.g:132:2: ( ( rule__PortElement__Alternatives ) )
            {
            // InternalElkGraphJson.g:132:2: ( ( rule__PortElement__Alternatives ) )
            // InternalElkGraphJson.g:133:3: ( rule__PortElement__Alternatives )
            {
             before(grammarAccess.getPortElementAccess().getAlternatives()); 
            // InternalElkGraphJson.g:134:3: ( rule__PortElement__Alternatives )
            // InternalElkGraphJson.g:134:4: rule__PortElement__Alternatives
            {
            pushFollow(FOLLOW_2);
            rule__PortElement__Alternatives();

            state._fsp--;


            }

             after(grammarAccess.getPortElementAccess().getAlternatives()); 

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
    // $ANTLR end "rulePortElement"


    // $ANTLR start "entryRuleElkLabel"
    // InternalElkGraphJson.g:143:1: entryRuleElkLabel : ruleElkLabel EOF ;
    public final void entryRuleElkLabel() throws RecognitionException {
        try {
            // InternalElkGraphJson.g:144:1: ( ruleElkLabel EOF )
            // InternalElkGraphJson.g:145:1: ruleElkLabel EOF
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
    // InternalElkGraphJson.g:152:1: ruleElkLabel : ( ( rule__ElkLabel__Group__0 ) ) ;
    public final void ruleElkLabel() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:156:2: ( ( ( rule__ElkLabel__Group__0 ) ) )
            // InternalElkGraphJson.g:157:2: ( ( rule__ElkLabel__Group__0 ) )
            {
            // InternalElkGraphJson.g:157:2: ( ( rule__ElkLabel__Group__0 ) )
            // InternalElkGraphJson.g:158:3: ( rule__ElkLabel__Group__0 )
            {
             before(grammarAccess.getElkLabelAccess().getGroup()); 
            // InternalElkGraphJson.g:159:3: ( rule__ElkLabel__Group__0 )
            // InternalElkGraphJson.g:159:4: rule__ElkLabel__Group__0
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


    // $ANTLR start "ruleLabelElement"
    // InternalElkGraphJson.g:169:1: ruleLabelElement : ( ( rule__LabelElement__Alternatives ) ) ;
    public final void ruleLabelElement() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:173:2: ( ( ( rule__LabelElement__Alternatives ) ) )
            // InternalElkGraphJson.g:174:2: ( ( rule__LabelElement__Alternatives ) )
            {
            // InternalElkGraphJson.g:174:2: ( ( rule__LabelElement__Alternatives ) )
            // InternalElkGraphJson.g:175:3: ( rule__LabelElement__Alternatives )
            {
             before(grammarAccess.getLabelElementAccess().getAlternatives()); 
            // InternalElkGraphJson.g:176:3: ( rule__LabelElement__Alternatives )
            // InternalElkGraphJson.g:176:4: rule__LabelElement__Alternatives
            {
            pushFollow(FOLLOW_2);
            rule__LabelElement__Alternatives();

            state._fsp--;


            }

             after(grammarAccess.getLabelElementAccess().getAlternatives()); 

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
    // $ANTLR end "ruleLabelElement"


    // $ANTLR start "entryRuleElkEdge"
    // InternalElkGraphJson.g:185:1: entryRuleElkEdge : ruleElkEdge EOF ;
    public final void entryRuleElkEdge() throws RecognitionException {
        try {
            // InternalElkGraphJson.g:186:1: ( ruleElkEdge EOF )
            // InternalElkGraphJson.g:187:1: ruleElkEdge EOF
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
    // InternalElkGraphJson.g:194:1: ruleElkEdge : ( ( rule__ElkEdge__Group__0 ) ) ;
    public final void ruleElkEdge() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:198:2: ( ( ( rule__ElkEdge__Group__0 ) ) )
            // InternalElkGraphJson.g:199:2: ( ( rule__ElkEdge__Group__0 ) )
            {
            // InternalElkGraphJson.g:199:2: ( ( rule__ElkEdge__Group__0 ) )
            // InternalElkGraphJson.g:200:3: ( rule__ElkEdge__Group__0 )
            {
             before(grammarAccess.getElkEdgeAccess().getGroup()); 
            // InternalElkGraphJson.g:201:3: ( rule__ElkEdge__Group__0 )
            // InternalElkGraphJson.g:201:4: rule__ElkEdge__Group__0
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


    // $ANTLR start "ruleEdgeElement"
    // InternalElkGraphJson.g:211:1: ruleEdgeElement : ( ( rule__EdgeElement__Alternatives ) ) ;
    public final void ruleEdgeElement() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:215:2: ( ( ( rule__EdgeElement__Alternatives ) ) )
            // InternalElkGraphJson.g:216:2: ( ( rule__EdgeElement__Alternatives ) )
            {
            // InternalElkGraphJson.g:216:2: ( ( rule__EdgeElement__Alternatives ) )
            // InternalElkGraphJson.g:217:3: ( rule__EdgeElement__Alternatives )
            {
             before(grammarAccess.getEdgeElementAccess().getAlternatives()); 
            // InternalElkGraphJson.g:218:3: ( rule__EdgeElement__Alternatives )
            // InternalElkGraphJson.g:218:4: rule__EdgeElement__Alternatives
            {
            pushFollow(FOLLOW_2);
            rule__EdgeElement__Alternatives();

            state._fsp--;


            }

             after(grammarAccess.getEdgeElementAccess().getAlternatives()); 

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
    // $ANTLR end "ruleEdgeElement"


    // $ANTLR start "ruleElkEdgeSources"
    // InternalElkGraphJson.g:228:1: ruleElkEdgeSources : ( ( rule__ElkEdgeSources__Group__0 ) ) ;
    public final void ruleElkEdgeSources() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:232:2: ( ( ( rule__ElkEdgeSources__Group__0 ) ) )
            // InternalElkGraphJson.g:233:2: ( ( rule__ElkEdgeSources__Group__0 ) )
            {
            // InternalElkGraphJson.g:233:2: ( ( rule__ElkEdgeSources__Group__0 ) )
            // InternalElkGraphJson.g:234:3: ( rule__ElkEdgeSources__Group__0 )
            {
             before(grammarAccess.getElkEdgeSourcesAccess().getGroup()); 
            // InternalElkGraphJson.g:235:3: ( rule__ElkEdgeSources__Group__0 )
            // InternalElkGraphJson.g:235:4: rule__ElkEdgeSources__Group__0
            {
            pushFollow(FOLLOW_2);
            rule__ElkEdgeSources__Group__0();

            state._fsp--;


            }

             after(grammarAccess.getElkEdgeSourcesAccess().getGroup()); 

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
    // $ANTLR end "ruleElkEdgeSources"


    // $ANTLR start "ruleElkEdgeTargets"
    // InternalElkGraphJson.g:245:1: ruleElkEdgeTargets : ( ( rule__ElkEdgeTargets__Group__0 ) ) ;
    public final void ruleElkEdgeTargets() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:249:2: ( ( ( rule__ElkEdgeTargets__Group__0 ) ) )
            // InternalElkGraphJson.g:250:2: ( ( rule__ElkEdgeTargets__Group__0 ) )
            {
            // InternalElkGraphJson.g:250:2: ( ( rule__ElkEdgeTargets__Group__0 ) )
            // InternalElkGraphJson.g:251:3: ( rule__ElkEdgeTargets__Group__0 )
            {
             before(grammarAccess.getElkEdgeTargetsAccess().getGroup()); 
            // InternalElkGraphJson.g:252:3: ( rule__ElkEdgeTargets__Group__0 )
            // InternalElkGraphJson.g:252:4: rule__ElkEdgeTargets__Group__0
            {
            pushFollow(FOLLOW_2);
            rule__ElkEdgeTargets__Group__0();

            state._fsp--;


            }

             after(grammarAccess.getElkEdgeTargetsAccess().getGroup()); 

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
    // $ANTLR end "ruleElkEdgeTargets"


    // $ANTLR start "ruleElkId"
    // InternalElkGraphJson.g:262:1: ruleElkId : ( ( rule__ElkId__Group__0 ) ) ;
    public final void ruleElkId() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:266:2: ( ( ( rule__ElkId__Group__0 ) ) )
            // InternalElkGraphJson.g:267:2: ( ( rule__ElkId__Group__0 ) )
            {
            // InternalElkGraphJson.g:267:2: ( ( rule__ElkId__Group__0 ) )
            // InternalElkGraphJson.g:268:3: ( rule__ElkId__Group__0 )
            {
             before(grammarAccess.getElkIdAccess().getGroup()); 
            // InternalElkGraphJson.g:269:3: ( rule__ElkId__Group__0 )
            // InternalElkGraphJson.g:269:4: rule__ElkId__Group__0
            {
            pushFollow(FOLLOW_2);
            rule__ElkId__Group__0();

            state._fsp--;


            }

             after(grammarAccess.getElkIdAccess().getGroup()); 

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
    // $ANTLR end "ruleElkId"


    // $ANTLR start "ruleElkNodeChildren"
    // InternalElkGraphJson.g:279:1: ruleElkNodeChildren : ( ( rule__ElkNodeChildren__Group__0 ) ) ;
    public final void ruleElkNodeChildren() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:283:2: ( ( ( rule__ElkNodeChildren__Group__0 ) ) )
            // InternalElkGraphJson.g:284:2: ( ( rule__ElkNodeChildren__Group__0 ) )
            {
            // InternalElkGraphJson.g:284:2: ( ( rule__ElkNodeChildren__Group__0 ) )
            // InternalElkGraphJson.g:285:3: ( rule__ElkNodeChildren__Group__0 )
            {
             before(grammarAccess.getElkNodeChildrenAccess().getGroup()); 
            // InternalElkGraphJson.g:286:3: ( rule__ElkNodeChildren__Group__0 )
            // InternalElkGraphJson.g:286:4: rule__ElkNodeChildren__Group__0
            {
            pushFollow(FOLLOW_2);
            rule__ElkNodeChildren__Group__0();

            state._fsp--;


            }

             after(grammarAccess.getElkNodeChildrenAccess().getGroup()); 

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
    // $ANTLR end "ruleElkNodeChildren"


    // $ANTLR start "ruleElkNodePorts"
    // InternalElkGraphJson.g:296:1: ruleElkNodePorts : ( ( rule__ElkNodePorts__Group__0 ) ) ;
    public final void ruleElkNodePorts() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:300:2: ( ( ( rule__ElkNodePorts__Group__0 ) ) )
            // InternalElkGraphJson.g:301:2: ( ( rule__ElkNodePorts__Group__0 ) )
            {
            // InternalElkGraphJson.g:301:2: ( ( rule__ElkNodePorts__Group__0 ) )
            // InternalElkGraphJson.g:302:3: ( rule__ElkNodePorts__Group__0 )
            {
             before(grammarAccess.getElkNodePortsAccess().getGroup()); 
            // InternalElkGraphJson.g:303:3: ( rule__ElkNodePorts__Group__0 )
            // InternalElkGraphJson.g:303:4: rule__ElkNodePorts__Group__0
            {
            pushFollow(FOLLOW_2);
            rule__ElkNodePorts__Group__0();

            state._fsp--;


            }

             after(grammarAccess.getElkNodePortsAccess().getGroup()); 

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
    // $ANTLR end "ruleElkNodePorts"


    // $ANTLR start "ruleElkNodeEdges"
    // InternalElkGraphJson.g:313:1: ruleElkNodeEdges : ( ( rule__ElkNodeEdges__Group__0 ) ) ;
    public final void ruleElkNodeEdges() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:317:2: ( ( ( rule__ElkNodeEdges__Group__0 ) ) )
            // InternalElkGraphJson.g:318:2: ( ( rule__ElkNodeEdges__Group__0 ) )
            {
            // InternalElkGraphJson.g:318:2: ( ( rule__ElkNodeEdges__Group__0 ) )
            // InternalElkGraphJson.g:319:3: ( rule__ElkNodeEdges__Group__0 )
            {
             before(grammarAccess.getElkNodeEdgesAccess().getGroup()); 
            // InternalElkGraphJson.g:320:3: ( rule__ElkNodeEdges__Group__0 )
            // InternalElkGraphJson.g:320:4: rule__ElkNodeEdges__Group__0
            {
            pushFollow(FOLLOW_2);
            rule__ElkNodeEdges__Group__0();

            state._fsp--;


            }

             after(grammarAccess.getElkNodeEdgesAccess().getGroup()); 

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
    // $ANTLR end "ruleElkNodeEdges"


    // $ANTLR start "ruleElkGraphElementLabels"
    // InternalElkGraphJson.g:330:1: ruleElkGraphElementLabels : ( ( rule__ElkGraphElementLabels__Group__0 ) ) ;
    public final void ruleElkGraphElementLabels() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:334:2: ( ( ( rule__ElkGraphElementLabels__Group__0 ) ) )
            // InternalElkGraphJson.g:335:2: ( ( rule__ElkGraphElementLabels__Group__0 ) )
            {
            // InternalElkGraphJson.g:335:2: ( ( rule__ElkGraphElementLabels__Group__0 ) )
            // InternalElkGraphJson.g:336:3: ( rule__ElkGraphElementLabels__Group__0 )
            {
             before(grammarAccess.getElkGraphElementLabelsAccess().getGroup()); 
            // InternalElkGraphJson.g:337:3: ( rule__ElkGraphElementLabels__Group__0 )
            // InternalElkGraphJson.g:337:4: rule__ElkGraphElementLabels__Group__0
            {
            pushFollow(FOLLOW_2);
            rule__ElkGraphElementLabels__Group__0();

            state._fsp--;


            }

             after(grammarAccess.getElkGraphElementLabelsAccess().getGroup()); 

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
    // $ANTLR end "ruleElkGraphElementLabels"


    // $ANTLR start "ruleElkGraphElementProperties"
    // InternalElkGraphJson.g:347:1: ruleElkGraphElementProperties : ( ( rule__ElkGraphElementProperties__Group__0 ) ) ;
    public final void ruleElkGraphElementProperties() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:351:2: ( ( ( rule__ElkGraphElementProperties__Group__0 ) ) )
            // InternalElkGraphJson.g:352:2: ( ( rule__ElkGraphElementProperties__Group__0 ) )
            {
            // InternalElkGraphJson.g:352:2: ( ( rule__ElkGraphElementProperties__Group__0 ) )
            // InternalElkGraphJson.g:353:3: ( rule__ElkGraphElementProperties__Group__0 )
            {
             before(grammarAccess.getElkGraphElementPropertiesAccess().getGroup()); 
            // InternalElkGraphJson.g:354:3: ( rule__ElkGraphElementProperties__Group__0 )
            // InternalElkGraphJson.g:354:4: rule__ElkGraphElementProperties__Group__0
            {
            pushFollow(FOLLOW_2);
            rule__ElkGraphElementProperties__Group__0();

            state._fsp--;


            }

             after(grammarAccess.getElkGraphElementPropertiesAccess().getGroup()); 

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
    // $ANTLR end "ruleElkGraphElementProperties"


    // $ANTLR start "ruleShapeElement"
    // InternalElkGraphJson.g:364:1: ruleShapeElement : ( ( rule__ShapeElement__Alternatives ) ) ;
    public final void ruleShapeElement() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:368:2: ( ( ( rule__ShapeElement__Alternatives ) ) )
            // InternalElkGraphJson.g:369:2: ( ( rule__ShapeElement__Alternatives ) )
            {
            // InternalElkGraphJson.g:369:2: ( ( rule__ShapeElement__Alternatives ) )
            // InternalElkGraphJson.g:370:3: ( rule__ShapeElement__Alternatives )
            {
             before(grammarAccess.getShapeElementAccess().getAlternatives()); 
            // InternalElkGraphJson.g:371:3: ( rule__ShapeElement__Alternatives )
            // InternalElkGraphJson.g:371:4: rule__ShapeElement__Alternatives
            {
            pushFollow(FOLLOW_2);
            rule__ShapeElement__Alternatives();

            state._fsp--;


            }

             after(grammarAccess.getShapeElementAccess().getAlternatives()); 

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
    // $ANTLR end "ruleShapeElement"


    // $ANTLR start "entryRuleProperty"
    // InternalElkGraphJson.g:380:1: entryRuleProperty : ruleProperty EOF ;
    public final void entryRuleProperty() throws RecognitionException {
        try {
            // InternalElkGraphJson.g:381:1: ( ruleProperty EOF )
            // InternalElkGraphJson.g:382:1: ruleProperty EOF
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
    // InternalElkGraphJson.g:389:1: ruleProperty : ( ( rule__Property__Group__0 ) ) ;
    public final void ruleProperty() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:393:2: ( ( ( rule__Property__Group__0 ) ) )
            // InternalElkGraphJson.g:394:2: ( ( rule__Property__Group__0 ) )
            {
            // InternalElkGraphJson.g:394:2: ( ( rule__Property__Group__0 ) )
            // InternalElkGraphJson.g:395:3: ( rule__Property__Group__0 )
            {
             before(grammarAccess.getPropertyAccess().getGroup()); 
            // InternalElkGraphJson.g:396:3: ( rule__Property__Group__0 )
            // InternalElkGraphJson.g:396:4: rule__Property__Group__0
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
    // InternalElkGraphJson.g:405:1: entryRulePropertyKey : rulePropertyKey EOF ;
    public final void entryRulePropertyKey() throws RecognitionException {
         
        	HiddenTokens myHiddenTokenState = ((XtextTokenStream)input).setHiddenTokens();

        try {
            // InternalElkGraphJson.g:409:1: ( rulePropertyKey EOF )
            // InternalElkGraphJson.g:410:1: rulePropertyKey EOF
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
    // InternalElkGraphJson.g:420:1: rulePropertyKey : ( ( rule__PropertyKey__Alternatives ) ) ;
    public final void rulePropertyKey() throws RecognitionException {

        		HiddenTokens myHiddenTokenState = ((XtextTokenStream)input).setHiddenTokens();
        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:425:2: ( ( ( rule__PropertyKey__Alternatives ) ) )
            // InternalElkGraphJson.g:426:2: ( ( rule__PropertyKey__Alternatives ) )
            {
            // InternalElkGraphJson.g:426:2: ( ( rule__PropertyKey__Alternatives ) )
            // InternalElkGraphJson.g:427:3: ( rule__PropertyKey__Alternatives )
            {
             before(grammarAccess.getPropertyKeyAccess().getAlternatives()); 
            // InternalElkGraphJson.g:428:3: ( rule__PropertyKey__Alternatives )
            // InternalElkGraphJson.g:428:4: rule__PropertyKey__Alternatives
            {
            pushFollow(FOLLOW_2);
            rule__PropertyKey__Alternatives();

            state._fsp--;


            }

             after(grammarAccess.getPropertyKeyAccess().getAlternatives()); 

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
    // InternalElkGraphJson.g:438:1: entryRuleStringValue : ruleStringValue EOF ;
    public final void entryRuleStringValue() throws RecognitionException {
        try {
            // InternalElkGraphJson.g:439:1: ( ruleStringValue EOF )
            // InternalElkGraphJson.g:440:1: ruleStringValue EOF
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
    // InternalElkGraphJson.g:447:1: ruleStringValue : ( RULE_STRING ) ;
    public final void ruleStringValue() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:451:2: ( ( RULE_STRING ) )
            // InternalElkGraphJson.g:452:2: ( RULE_STRING )
            {
            // InternalElkGraphJson.g:452:2: ( RULE_STRING )
            // InternalElkGraphJson.g:453:3: RULE_STRING
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


    // $ANTLR start "entryRuleNumberValue"
    // InternalElkGraphJson.g:463:1: entryRuleNumberValue : ruleNumberValue EOF ;
    public final void entryRuleNumberValue() throws RecognitionException {
        try {
            // InternalElkGraphJson.g:464:1: ( ruleNumberValue EOF )
            // InternalElkGraphJson.g:465:1: ruleNumberValue EOF
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
    // InternalElkGraphJson.g:472:1: ruleNumberValue : ( ( rule__NumberValue__Alternatives ) ) ;
    public final void ruleNumberValue() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:476:2: ( ( ( rule__NumberValue__Alternatives ) ) )
            // InternalElkGraphJson.g:477:2: ( ( rule__NumberValue__Alternatives ) )
            {
            // InternalElkGraphJson.g:477:2: ( ( rule__NumberValue__Alternatives ) )
            // InternalElkGraphJson.g:478:3: ( rule__NumberValue__Alternatives )
            {
             before(grammarAccess.getNumberValueAccess().getAlternatives()); 
            // InternalElkGraphJson.g:479:3: ( rule__NumberValue__Alternatives )
            // InternalElkGraphJson.g:479:4: rule__NumberValue__Alternatives
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
    // InternalElkGraphJson.g:488:1: entryRuleBooleanValue : ruleBooleanValue EOF ;
    public final void entryRuleBooleanValue() throws RecognitionException {
        try {
            // InternalElkGraphJson.g:489:1: ( ruleBooleanValue EOF )
            // InternalElkGraphJson.g:490:1: ruleBooleanValue EOF
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
    // InternalElkGraphJson.g:497:1: ruleBooleanValue : ( ( rule__BooleanValue__Alternatives ) ) ;
    public final void ruleBooleanValue() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:501:2: ( ( ( rule__BooleanValue__Alternatives ) ) )
            // InternalElkGraphJson.g:502:2: ( ( rule__BooleanValue__Alternatives ) )
            {
            // InternalElkGraphJson.g:502:2: ( ( rule__BooleanValue__Alternatives ) )
            // InternalElkGraphJson.g:503:3: ( rule__BooleanValue__Alternatives )
            {
             before(grammarAccess.getBooleanValueAccess().getAlternatives()); 
            // InternalElkGraphJson.g:504:3: ( rule__BooleanValue__Alternatives )
            // InternalElkGraphJson.g:504:4: rule__BooleanValue__Alternatives
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


    // $ANTLR start "entryRuleNumber"
    // InternalElkGraphJson.g:513:1: entryRuleNumber : ruleNumber EOF ;
    public final void entryRuleNumber() throws RecognitionException {
        try {
            // InternalElkGraphJson.g:514:1: ( ruleNumber EOF )
            // InternalElkGraphJson.g:515:1: ruleNumber EOF
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
    // InternalElkGraphJson.g:522:1: ruleNumber : ( ( rule__Number__Alternatives ) ) ;
    public final void ruleNumber() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:526:2: ( ( ( rule__Number__Alternatives ) ) )
            // InternalElkGraphJson.g:527:2: ( ( rule__Number__Alternatives ) )
            {
            // InternalElkGraphJson.g:527:2: ( ( rule__Number__Alternatives ) )
            // InternalElkGraphJson.g:528:3: ( rule__Number__Alternatives )
            {
             before(grammarAccess.getNumberAccess().getAlternatives()); 
            // InternalElkGraphJson.g:529:3: ( rule__Number__Alternatives )
            // InternalElkGraphJson.g:529:4: rule__Number__Alternatives
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


    // $ANTLR start "entryRuleJsonObject"
    // InternalElkGraphJson.g:538:1: entryRuleJsonObject : ruleJsonObject EOF ;
    public final void entryRuleJsonObject() throws RecognitionException {
        try {
            // InternalElkGraphJson.g:539:1: ( ruleJsonObject EOF )
            // InternalElkGraphJson.g:540:1: ruleJsonObject EOF
            {
             before(grammarAccess.getJsonObjectRule()); 
            pushFollow(FOLLOW_1);
            ruleJsonObject();

            state._fsp--;

             after(grammarAccess.getJsonObjectRule()); 
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
    // $ANTLR end "entryRuleJsonObject"


    // $ANTLR start "ruleJsonObject"
    // InternalElkGraphJson.g:547:1: ruleJsonObject : ( ( rule__JsonObject__Group__0 ) ) ;
    public final void ruleJsonObject() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:551:2: ( ( ( rule__JsonObject__Group__0 ) ) )
            // InternalElkGraphJson.g:552:2: ( ( rule__JsonObject__Group__0 ) )
            {
            // InternalElkGraphJson.g:552:2: ( ( rule__JsonObject__Group__0 ) )
            // InternalElkGraphJson.g:553:3: ( rule__JsonObject__Group__0 )
            {
             before(grammarAccess.getJsonObjectAccess().getGroup()); 
            // InternalElkGraphJson.g:554:3: ( rule__JsonObject__Group__0 )
            // InternalElkGraphJson.g:554:4: rule__JsonObject__Group__0
            {
            pushFollow(FOLLOW_2);
            rule__JsonObject__Group__0();

            state._fsp--;


            }

             after(grammarAccess.getJsonObjectAccess().getGroup()); 

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
    // $ANTLR end "ruleJsonObject"


    // $ANTLR start "entryRuleJsonArray"
    // InternalElkGraphJson.g:563:1: entryRuleJsonArray : ruleJsonArray EOF ;
    public final void entryRuleJsonArray() throws RecognitionException {
        try {
            // InternalElkGraphJson.g:564:1: ( ruleJsonArray EOF )
            // InternalElkGraphJson.g:565:1: ruleJsonArray EOF
            {
             before(grammarAccess.getJsonArrayRule()); 
            pushFollow(FOLLOW_1);
            ruleJsonArray();

            state._fsp--;

             after(grammarAccess.getJsonArrayRule()); 
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
    // $ANTLR end "entryRuleJsonArray"


    // $ANTLR start "ruleJsonArray"
    // InternalElkGraphJson.g:572:1: ruleJsonArray : ( ( rule__JsonArray__Group__0 ) ) ;
    public final void ruleJsonArray() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:576:2: ( ( ( rule__JsonArray__Group__0 ) ) )
            // InternalElkGraphJson.g:577:2: ( ( rule__JsonArray__Group__0 ) )
            {
            // InternalElkGraphJson.g:577:2: ( ( rule__JsonArray__Group__0 ) )
            // InternalElkGraphJson.g:578:3: ( rule__JsonArray__Group__0 )
            {
             before(grammarAccess.getJsonArrayAccess().getGroup()); 
            // InternalElkGraphJson.g:579:3: ( rule__JsonArray__Group__0 )
            // InternalElkGraphJson.g:579:4: rule__JsonArray__Group__0
            {
            pushFollow(FOLLOW_2);
            rule__JsonArray__Group__0();

            state._fsp--;


            }

             after(grammarAccess.getJsonArrayAccess().getGroup()); 

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
    // $ANTLR end "ruleJsonArray"


    // $ANTLR start "entryRuleJsonMember"
    // InternalElkGraphJson.g:588:1: entryRuleJsonMember : ruleJsonMember EOF ;
    public final void entryRuleJsonMember() throws RecognitionException {
        try {
            // InternalElkGraphJson.g:589:1: ( ruleJsonMember EOF )
            // InternalElkGraphJson.g:590:1: ruleJsonMember EOF
            {
             before(grammarAccess.getJsonMemberRule()); 
            pushFollow(FOLLOW_1);
            ruleJsonMember();

            state._fsp--;

             after(grammarAccess.getJsonMemberRule()); 
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
    // $ANTLR end "entryRuleJsonMember"


    // $ANTLR start "ruleJsonMember"
    // InternalElkGraphJson.g:597:1: ruleJsonMember : ( ( rule__JsonMember__Group__0 ) ) ;
    public final void ruleJsonMember() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:601:2: ( ( ( rule__JsonMember__Group__0 ) ) )
            // InternalElkGraphJson.g:602:2: ( ( rule__JsonMember__Group__0 ) )
            {
            // InternalElkGraphJson.g:602:2: ( ( rule__JsonMember__Group__0 ) )
            // InternalElkGraphJson.g:603:3: ( rule__JsonMember__Group__0 )
            {
             before(grammarAccess.getJsonMemberAccess().getGroup()); 
            // InternalElkGraphJson.g:604:3: ( rule__JsonMember__Group__0 )
            // InternalElkGraphJson.g:604:4: rule__JsonMember__Group__0
            {
            pushFollow(FOLLOW_2);
            rule__JsonMember__Group__0();

            state._fsp--;


            }

             after(grammarAccess.getJsonMemberAccess().getGroup()); 

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
    // $ANTLR end "ruleJsonMember"


    // $ANTLR start "entryRuleJsonValue"
    // InternalElkGraphJson.g:613:1: entryRuleJsonValue : ruleJsonValue EOF ;
    public final void entryRuleJsonValue() throws RecognitionException {
        try {
            // InternalElkGraphJson.g:614:1: ( ruleJsonValue EOF )
            // InternalElkGraphJson.g:615:1: ruleJsonValue EOF
            {
             before(grammarAccess.getJsonValueRule()); 
            pushFollow(FOLLOW_1);
            ruleJsonValue();

            state._fsp--;

             after(grammarAccess.getJsonValueRule()); 
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
    // $ANTLR end "entryRuleJsonValue"


    // $ANTLR start "ruleJsonValue"
    // InternalElkGraphJson.g:622:1: ruleJsonValue : ( ( rule__JsonValue__Alternatives ) ) ;
    public final void ruleJsonValue() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:626:2: ( ( ( rule__JsonValue__Alternatives ) ) )
            // InternalElkGraphJson.g:627:2: ( ( rule__JsonValue__Alternatives ) )
            {
            // InternalElkGraphJson.g:627:2: ( ( rule__JsonValue__Alternatives ) )
            // InternalElkGraphJson.g:628:3: ( rule__JsonValue__Alternatives )
            {
             before(grammarAccess.getJsonValueAccess().getAlternatives()); 
            // InternalElkGraphJson.g:629:3: ( rule__JsonValue__Alternatives )
            // InternalElkGraphJson.g:629:4: rule__JsonValue__Alternatives
            {
            pushFollow(FOLLOW_2);
            rule__JsonValue__Alternatives();

            state._fsp--;


            }

             after(grammarAccess.getJsonValueAccess().getAlternatives()); 

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
    // $ANTLR end "ruleJsonValue"


    // $ANTLR start "entryRuleKeyChildren"
    // InternalElkGraphJson.g:638:1: entryRuleKeyChildren : ruleKeyChildren EOF ;
    public final void entryRuleKeyChildren() throws RecognitionException {
        try {
            // InternalElkGraphJson.g:639:1: ( ruleKeyChildren EOF )
            // InternalElkGraphJson.g:640:1: ruleKeyChildren EOF
            {
             before(grammarAccess.getKeyChildrenRule()); 
            pushFollow(FOLLOW_1);
            ruleKeyChildren();

            state._fsp--;

             after(grammarAccess.getKeyChildrenRule()); 
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
    // $ANTLR end "entryRuleKeyChildren"


    // $ANTLR start "ruleKeyChildren"
    // InternalElkGraphJson.g:647:1: ruleKeyChildren : ( ( rule__KeyChildren__Alternatives ) ) ;
    public final void ruleKeyChildren() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:651:2: ( ( ( rule__KeyChildren__Alternatives ) ) )
            // InternalElkGraphJson.g:652:2: ( ( rule__KeyChildren__Alternatives ) )
            {
            // InternalElkGraphJson.g:652:2: ( ( rule__KeyChildren__Alternatives ) )
            // InternalElkGraphJson.g:653:3: ( rule__KeyChildren__Alternatives )
            {
             before(grammarAccess.getKeyChildrenAccess().getAlternatives()); 
            // InternalElkGraphJson.g:654:3: ( rule__KeyChildren__Alternatives )
            // InternalElkGraphJson.g:654:4: rule__KeyChildren__Alternatives
            {
            pushFollow(FOLLOW_2);
            rule__KeyChildren__Alternatives();

            state._fsp--;


            }

             after(grammarAccess.getKeyChildrenAccess().getAlternatives()); 

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
    // $ANTLR end "ruleKeyChildren"


    // $ANTLR start "entryRuleKeyPorts"
    // InternalElkGraphJson.g:663:1: entryRuleKeyPorts : ruleKeyPorts EOF ;
    public final void entryRuleKeyPorts() throws RecognitionException {
        try {
            // InternalElkGraphJson.g:664:1: ( ruleKeyPorts EOF )
            // InternalElkGraphJson.g:665:1: ruleKeyPorts EOF
            {
             before(grammarAccess.getKeyPortsRule()); 
            pushFollow(FOLLOW_1);
            ruleKeyPorts();

            state._fsp--;

             after(grammarAccess.getKeyPortsRule()); 
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
    // $ANTLR end "entryRuleKeyPorts"


    // $ANTLR start "ruleKeyPorts"
    // InternalElkGraphJson.g:672:1: ruleKeyPorts : ( ( rule__KeyPorts__Alternatives ) ) ;
    public final void ruleKeyPorts() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:676:2: ( ( ( rule__KeyPorts__Alternatives ) ) )
            // InternalElkGraphJson.g:677:2: ( ( rule__KeyPorts__Alternatives ) )
            {
            // InternalElkGraphJson.g:677:2: ( ( rule__KeyPorts__Alternatives ) )
            // InternalElkGraphJson.g:678:3: ( rule__KeyPorts__Alternatives )
            {
             before(grammarAccess.getKeyPortsAccess().getAlternatives()); 
            // InternalElkGraphJson.g:679:3: ( rule__KeyPorts__Alternatives )
            // InternalElkGraphJson.g:679:4: rule__KeyPorts__Alternatives
            {
            pushFollow(FOLLOW_2);
            rule__KeyPorts__Alternatives();

            state._fsp--;


            }

             after(grammarAccess.getKeyPortsAccess().getAlternatives()); 

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
    // $ANTLR end "ruleKeyPorts"


    // $ANTLR start "entryRuleKeyLabels"
    // InternalElkGraphJson.g:688:1: entryRuleKeyLabels : ruleKeyLabels EOF ;
    public final void entryRuleKeyLabels() throws RecognitionException {
        try {
            // InternalElkGraphJson.g:689:1: ( ruleKeyLabels EOF )
            // InternalElkGraphJson.g:690:1: ruleKeyLabels EOF
            {
             before(grammarAccess.getKeyLabelsRule()); 
            pushFollow(FOLLOW_1);
            ruleKeyLabels();

            state._fsp--;

             after(grammarAccess.getKeyLabelsRule()); 
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
    // $ANTLR end "entryRuleKeyLabels"


    // $ANTLR start "ruleKeyLabels"
    // InternalElkGraphJson.g:697:1: ruleKeyLabels : ( ( rule__KeyLabels__Alternatives ) ) ;
    public final void ruleKeyLabels() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:701:2: ( ( ( rule__KeyLabels__Alternatives ) ) )
            // InternalElkGraphJson.g:702:2: ( ( rule__KeyLabels__Alternatives ) )
            {
            // InternalElkGraphJson.g:702:2: ( ( rule__KeyLabels__Alternatives ) )
            // InternalElkGraphJson.g:703:3: ( rule__KeyLabels__Alternatives )
            {
             before(grammarAccess.getKeyLabelsAccess().getAlternatives()); 
            // InternalElkGraphJson.g:704:3: ( rule__KeyLabels__Alternatives )
            // InternalElkGraphJson.g:704:4: rule__KeyLabels__Alternatives
            {
            pushFollow(FOLLOW_2);
            rule__KeyLabels__Alternatives();

            state._fsp--;


            }

             after(grammarAccess.getKeyLabelsAccess().getAlternatives()); 

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
    // $ANTLR end "ruleKeyLabels"


    // $ANTLR start "entryRuleKeyEdges"
    // InternalElkGraphJson.g:713:1: entryRuleKeyEdges : ruleKeyEdges EOF ;
    public final void entryRuleKeyEdges() throws RecognitionException {
        try {
            // InternalElkGraphJson.g:714:1: ( ruleKeyEdges EOF )
            // InternalElkGraphJson.g:715:1: ruleKeyEdges EOF
            {
             before(grammarAccess.getKeyEdgesRule()); 
            pushFollow(FOLLOW_1);
            ruleKeyEdges();

            state._fsp--;

             after(grammarAccess.getKeyEdgesRule()); 
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
    // $ANTLR end "entryRuleKeyEdges"


    // $ANTLR start "ruleKeyEdges"
    // InternalElkGraphJson.g:722:1: ruleKeyEdges : ( ( rule__KeyEdges__Alternatives ) ) ;
    public final void ruleKeyEdges() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:726:2: ( ( ( rule__KeyEdges__Alternatives ) ) )
            // InternalElkGraphJson.g:727:2: ( ( rule__KeyEdges__Alternatives ) )
            {
            // InternalElkGraphJson.g:727:2: ( ( rule__KeyEdges__Alternatives ) )
            // InternalElkGraphJson.g:728:3: ( rule__KeyEdges__Alternatives )
            {
             before(grammarAccess.getKeyEdgesAccess().getAlternatives()); 
            // InternalElkGraphJson.g:729:3: ( rule__KeyEdges__Alternatives )
            // InternalElkGraphJson.g:729:4: rule__KeyEdges__Alternatives
            {
            pushFollow(FOLLOW_2);
            rule__KeyEdges__Alternatives();

            state._fsp--;


            }

             after(grammarAccess.getKeyEdgesAccess().getAlternatives()); 

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
    // $ANTLR end "ruleKeyEdges"


    // $ANTLR start "entryRuleKeyLayoutOptions"
    // InternalElkGraphJson.g:738:1: entryRuleKeyLayoutOptions : ruleKeyLayoutOptions EOF ;
    public final void entryRuleKeyLayoutOptions() throws RecognitionException {
        try {
            // InternalElkGraphJson.g:739:1: ( ruleKeyLayoutOptions EOF )
            // InternalElkGraphJson.g:740:1: ruleKeyLayoutOptions EOF
            {
             before(grammarAccess.getKeyLayoutOptionsRule()); 
            pushFollow(FOLLOW_1);
            ruleKeyLayoutOptions();

            state._fsp--;

             after(grammarAccess.getKeyLayoutOptionsRule()); 
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
    // $ANTLR end "entryRuleKeyLayoutOptions"


    // $ANTLR start "ruleKeyLayoutOptions"
    // InternalElkGraphJson.g:747:1: ruleKeyLayoutOptions : ( ( rule__KeyLayoutOptions__Alternatives ) ) ;
    public final void ruleKeyLayoutOptions() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:751:2: ( ( ( rule__KeyLayoutOptions__Alternatives ) ) )
            // InternalElkGraphJson.g:752:2: ( ( rule__KeyLayoutOptions__Alternatives ) )
            {
            // InternalElkGraphJson.g:752:2: ( ( rule__KeyLayoutOptions__Alternatives ) )
            // InternalElkGraphJson.g:753:3: ( rule__KeyLayoutOptions__Alternatives )
            {
             before(grammarAccess.getKeyLayoutOptionsAccess().getAlternatives()); 
            // InternalElkGraphJson.g:754:3: ( rule__KeyLayoutOptions__Alternatives )
            // InternalElkGraphJson.g:754:4: rule__KeyLayoutOptions__Alternatives
            {
            pushFollow(FOLLOW_2);
            rule__KeyLayoutOptions__Alternatives();

            state._fsp--;


            }

             after(grammarAccess.getKeyLayoutOptionsAccess().getAlternatives()); 

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
    // $ANTLR end "ruleKeyLayoutOptions"


    // $ANTLR start "entryRuleKeyId"
    // InternalElkGraphJson.g:763:1: entryRuleKeyId : ruleKeyId EOF ;
    public final void entryRuleKeyId() throws RecognitionException {
        try {
            // InternalElkGraphJson.g:764:1: ( ruleKeyId EOF )
            // InternalElkGraphJson.g:765:1: ruleKeyId EOF
            {
             before(grammarAccess.getKeyIdRule()); 
            pushFollow(FOLLOW_1);
            ruleKeyId();

            state._fsp--;

             after(grammarAccess.getKeyIdRule()); 
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
    // $ANTLR end "entryRuleKeyId"


    // $ANTLR start "ruleKeyId"
    // InternalElkGraphJson.g:772:1: ruleKeyId : ( ( rule__KeyId__Alternatives ) ) ;
    public final void ruleKeyId() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:776:2: ( ( ( rule__KeyId__Alternatives ) ) )
            // InternalElkGraphJson.g:777:2: ( ( rule__KeyId__Alternatives ) )
            {
            // InternalElkGraphJson.g:777:2: ( ( rule__KeyId__Alternatives ) )
            // InternalElkGraphJson.g:778:3: ( rule__KeyId__Alternatives )
            {
             before(grammarAccess.getKeyIdAccess().getAlternatives()); 
            // InternalElkGraphJson.g:779:3: ( rule__KeyId__Alternatives )
            // InternalElkGraphJson.g:779:4: rule__KeyId__Alternatives
            {
            pushFollow(FOLLOW_2);
            rule__KeyId__Alternatives();

            state._fsp--;


            }

             after(grammarAccess.getKeyIdAccess().getAlternatives()); 

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
    // $ANTLR end "ruleKeyId"


    // $ANTLR start "entryRuleKeyX"
    // InternalElkGraphJson.g:788:1: entryRuleKeyX : ruleKeyX EOF ;
    public final void entryRuleKeyX() throws RecognitionException {
        try {
            // InternalElkGraphJson.g:789:1: ( ruleKeyX EOF )
            // InternalElkGraphJson.g:790:1: ruleKeyX EOF
            {
             before(grammarAccess.getKeyXRule()); 
            pushFollow(FOLLOW_1);
            ruleKeyX();

            state._fsp--;

             after(grammarAccess.getKeyXRule()); 
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
    // $ANTLR end "entryRuleKeyX"


    // $ANTLR start "ruleKeyX"
    // InternalElkGraphJson.g:797:1: ruleKeyX : ( ( rule__KeyX__Alternatives ) ) ;
    public final void ruleKeyX() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:801:2: ( ( ( rule__KeyX__Alternatives ) ) )
            // InternalElkGraphJson.g:802:2: ( ( rule__KeyX__Alternatives ) )
            {
            // InternalElkGraphJson.g:802:2: ( ( rule__KeyX__Alternatives ) )
            // InternalElkGraphJson.g:803:3: ( rule__KeyX__Alternatives )
            {
             before(grammarAccess.getKeyXAccess().getAlternatives()); 
            // InternalElkGraphJson.g:804:3: ( rule__KeyX__Alternatives )
            // InternalElkGraphJson.g:804:4: rule__KeyX__Alternatives
            {
            pushFollow(FOLLOW_2);
            rule__KeyX__Alternatives();

            state._fsp--;


            }

             after(grammarAccess.getKeyXAccess().getAlternatives()); 

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
    // $ANTLR end "ruleKeyX"


    // $ANTLR start "entryRuleKeyY"
    // InternalElkGraphJson.g:813:1: entryRuleKeyY : ruleKeyY EOF ;
    public final void entryRuleKeyY() throws RecognitionException {
        try {
            // InternalElkGraphJson.g:814:1: ( ruleKeyY EOF )
            // InternalElkGraphJson.g:815:1: ruleKeyY EOF
            {
             before(grammarAccess.getKeyYRule()); 
            pushFollow(FOLLOW_1);
            ruleKeyY();

            state._fsp--;

             after(grammarAccess.getKeyYRule()); 
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
    // $ANTLR end "entryRuleKeyY"


    // $ANTLR start "ruleKeyY"
    // InternalElkGraphJson.g:822:1: ruleKeyY : ( ( rule__KeyY__Alternatives ) ) ;
    public final void ruleKeyY() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:826:2: ( ( ( rule__KeyY__Alternatives ) ) )
            // InternalElkGraphJson.g:827:2: ( ( rule__KeyY__Alternatives ) )
            {
            // InternalElkGraphJson.g:827:2: ( ( rule__KeyY__Alternatives ) )
            // InternalElkGraphJson.g:828:3: ( rule__KeyY__Alternatives )
            {
             before(grammarAccess.getKeyYAccess().getAlternatives()); 
            // InternalElkGraphJson.g:829:3: ( rule__KeyY__Alternatives )
            // InternalElkGraphJson.g:829:4: rule__KeyY__Alternatives
            {
            pushFollow(FOLLOW_2);
            rule__KeyY__Alternatives();

            state._fsp--;


            }

             after(grammarAccess.getKeyYAccess().getAlternatives()); 

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
    // $ANTLR end "ruleKeyY"


    // $ANTLR start "entryRuleKeyWidth"
    // InternalElkGraphJson.g:838:1: entryRuleKeyWidth : ruleKeyWidth EOF ;
    public final void entryRuleKeyWidth() throws RecognitionException {
        try {
            // InternalElkGraphJson.g:839:1: ( ruleKeyWidth EOF )
            // InternalElkGraphJson.g:840:1: ruleKeyWidth EOF
            {
             before(grammarAccess.getKeyWidthRule()); 
            pushFollow(FOLLOW_1);
            ruleKeyWidth();

            state._fsp--;

             after(grammarAccess.getKeyWidthRule()); 
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
    // $ANTLR end "entryRuleKeyWidth"


    // $ANTLR start "ruleKeyWidth"
    // InternalElkGraphJson.g:847:1: ruleKeyWidth : ( ( rule__KeyWidth__Alternatives ) ) ;
    public final void ruleKeyWidth() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:851:2: ( ( ( rule__KeyWidth__Alternatives ) ) )
            // InternalElkGraphJson.g:852:2: ( ( rule__KeyWidth__Alternatives ) )
            {
            // InternalElkGraphJson.g:852:2: ( ( rule__KeyWidth__Alternatives ) )
            // InternalElkGraphJson.g:853:3: ( rule__KeyWidth__Alternatives )
            {
             before(grammarAccess.getKeyWidthAccess().getAlternatives()); 
            // InternalElkGraphJson.g:854:3: ( rule__KeyWidth__Alternatives )
            // InternalElkGraphJson.g:854:4: rule__KeyWidth__Alternatives
            {
            pushFollow(FOLLOW_2);
            rule__KeyWidth__Alternatives();

            state._fsp--;


            }

             after(grammarAccess.getKeyWidthAccess().getAlternatives()); 

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
    // $ANTLR end "ruleKeyWidth"


    // $ANTLR start "entryRuleKeyHeight"
    // InternalElkGraphJson.g:863:1: entryRuleKeyHeight : ruleKeyHeight EOF ;
    public final void entryRuleKeyHeight() throws RecognitionException {
        try {
            // InternalElkGraphJson.g:864:1: ( ruleKeyHeight EOF )
            // InternalElkGraphJson.g:865:1: ruleKeyHeight EOF
            {
             before(grammarAccess.getKeyHeightRule()); 
            pushFollow(FOLLOW_1);
            ruleKeyHeight();

            state._fsp--;

             after(grammarAccess.getKeyHeightRule()); 
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
    // $ANTLR end "entryRuleKeyHeight"


    // $ANTLR start "ruleKeyHeight"
    // InternalElkGraphJson.g:872:1: ruleKeyHeight : ( ( rule__KeyHeight__Alternatives ) ) ;
    public final void ruleKeyHeight() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:876:2: ( ( ( rule__KeyHeight__Alternatives ) ) )
            // InternalElkGraphJson.g:877:2: ( ( rule__KeyHeight__Alternatives ) )
            {
            // InternalElkGraphJson.g:877:2: ( ( rule__KeyHeight__Alternatives ) )
            // InternalElkGraphJson.g:878:3: ( rule__KeyHeight__Alternatives )
            {
             before(grammarAccess.getKeyHeightAccess().getAlternatives()); 
            // InternalElkGraphJson.g:879:3: ( rule__KeyHeight__Alternatives )
            // InternalElkGraphJson.g:879:4: rule__KeyHeight__Alternatives
            {
            pushFollow(FOLLOW_2);
            rule__KeyHeight__Alternatives();

            state._fsp--;


            }

             after(grammarAccess.getKeyHeightAccess().getAlternatives()); 

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
    // $ANTLR end "ruleKeyHeight"


    // $ANTLR start "entryRuleKeySources"
    // InternalElkGraphJson.g:888:1: entryRuleKeySources : ruleKeySources EOF ;
    public final void entryRuleKeySources() throws RecognitionException {
        try {
            // InternalElkGraphJson.g:889:1: ( ruleKeySources EOF )
            // InternalElkGraphJson.g:890:1: ruleKeySources EOF
            {
             before(grammarAccess.getKeySourcesRule()); 
            pushFollow(FOLLOW_1);
            ruleKeySources();

            state._fsp--;

             after(grammarAccess.getKeySourcesRule()); 
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
    // $ANTLR end "entryRuleKeySources"


    // $ANTLR start "ruleKeySources"
    // InternalElkGraphJson.g:897:1: ruleKeySources : ( ( rule__KeySources__Alternatives ) ) ;
    public final void ruleKeySources() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:901:2: ( ( ( rule__KeySources__Alternatives ) ) )
            // InternalElkGraphJson.g:902:2: ( ( rule__KeySources__Alternatives ) )
            {
            // InternalElkGraphJson.g:902:2: ( ( rule__KeySources__Alternatives ) )
            // InternalElkGraphJson.g:903:3: ( rule__KeySources__Alternatives )
            {
             before(grammarAccess.getKeySourcesAccess().getAlternatives()); 
            // InternalElkGraphJson.g:904:3: ( rule__KeySources__Alternatives )
            // InternalElkGraphJson.g:904:4: rule__KeySources__Alternatives
            {
            pushFollow(FOLLOW_2);
            rule__KeySources__Alternatives();

            state._fsp--;


            }

             after(grammarAccess.getKeySourcesAccess().getAlternatives()); 

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
    // $ANTLR end "ruleKeySources"


    // $ANTLR start "entryRuleKeyTargets"
    // InternalElkGraphJson.g:913:1: entryRuleKeyTargets : ruleKeyTargets EOF ;
    public final void entryRuleKeyTargets() throws RecognitionException {
        try {
            // InternalElkGraphJson.g:914:1: ( ruleKeyTargets EOF )
            // InternalElkGraphJson.g:915:1: ruleKeyTargets EOF
            {
             before(grammarAccess.getKeyTargetsRule()); 
            pushFollow(FOLLOW_1);
            ruleKeyTargets();

            state._fsp--;

             after(grammarAccess.getKeyTargetsRule()); 
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
    // $ANTLR end "entryRuleKeyTargets"


    // $ANTLR start "ruleKeyTargets"
    // InternalElkGraphJson.g:922:1: ruleKeyTargets : ( ( rule__KeyTargets__Alternatives ) ) ;
    public final void ruleKeyTargets() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:926:2: ( ( ( rule__KeyTargets__Alternatives ) ) )
            // InternalElkGraphJson.g:927:2: ( ( rule__KeyTargets__Alternatives ) )
            {
            // InternalElkGraphJson.g:927:2: ( ( rule__KeyTargets__Alternatives ) )
            // InternalElkGraphJson.g:928:3: ( rule__KeyTargets__Alternatives )
            {
             before(grammarAccess.getKeyTargetsAccess().getAlternatives()); 
            // InternalElkGraphJson.g:929:3: ( rule__KeyTargets__Alternatives )
            // InternalElkGraphJson.g:929:4: rule__KeyTargets__Alternatives
            {
            pushFollow(FOLLOW_2);
            rule__KeyTargets__Alternatives();

            state._fsp--;


            }

             after(grammarAccess.getKeyTargetsAccess().getAlternatives()); 

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
    // $ANTLR end "ruleKeyTargets"


    // $ANTLR start "entryRuleKeyText"
    // InternalElkGraphJson.g:938:1: entryRuleKeyText : ruleKeyText EOF ;
    public final void entryRuleKeyText() throws RecognitionException {
        try {
            // InternalElkGraphJson.g:939:1: ( ruleKeyText EOF )
            // InternalElkGraphJson.g:940:1: ruleKeyText EOF
            {
             before(grammarAccess.getKeyTextRule()); 
            pushFollow(FOLLOW_1);
            ruleKeyText();

            state._fsp--;

             after(grammarAccess.getKeyTextRule()); 
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
    // $ANTLR end "entryRuleKeyText"


    // $ANTLR start "ruleKeyText"
    // InternalElkGraphJson.g:947:1: ruleKeyText : ( ( rule__KeyText__Alternatives ) ) ;
    public final void ruleKeyText() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:951:2: ( ( ( rule__KeyText__Alternatives ) ) )
            // InternalElkGraphJson.g:952:2: ( ( rule__KeyText__Alternatives ) )
            {
            // InternalElkGraphJson.g:952:2: ( ( rule__KeyText__Alternatives ) )
            // InternalElkGraphJson.g:953:3: ( rule__KeyText__Alternatives )
            {
             before(grammarAccess.getKeyTextAccess().getAlternatives()); 
            // InternalElkGraphJson.g:954:3: ( rule__KeyText__Alternatives )
            // InternalElkGraphJson.g:954:4: rule__KeyText__Alternatives
            {
            pushFollow(FOLLOW_2);
            rule__KeyText__Alternatives();

            state._fsp--;


            }

             after(grammarAccess.getKeyTextAccess().getAlternatives()); 

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
    // $ANTLR end "ruleKeyText"


    // $ANTLR start "rule__NodeElement__Alternatives"
    // InternalElkGraphJson.g:962:1: rule__NodeElement__Alternatives : ( ( ruleElkId ) | ( ruleShapeElement ) | ( ( rule__NodeElement__Group_2__0 ) ) | ( ( rule__NodeElement__Group_3__0 ) ) | ( ( rule__NodeElement__Group_4__0 ) ) | ( ( rule__NodeElement__Group_5__0 ) ) | ( ( rule__NodeElement__Group_6__0 ) ) | ( ruleJsonMember ) );
    public final void rule__NodeElement__Alternatives() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:966:1: ( ( ruleElkId ) | ( ruleShapeElement ) | ( ( rule__NodeElement__Group_2__0 ) ) | ( ( rule__NodeElement__Group_3__0 ) ) | ( ( rule__NodeElement__Group_4__0 ) ) | ( ( rule__NodeElement__Group_5__0 ) ) | ( ( rule__NodeElement__Group_6__0 ) ) | ( ruleJsonMember ) )
            int alt1=8;
            switch ( input.LA(1) ) {
            case 34:
            case 35:
            case 36:
                {
                alt1=1;
                }
                break;
            case 37:
            case 38:
            case 39:
            case 40:
            case 41:
            case 42:
            case 43:
            case 44:
            case 45:
            case 46:
            case 47:
            case 48:
                {
                alt1=2;
                }
                break;
            case 16:
            case 17:
            case 18:
                {
                alt1=3;
                }
                break;
            case 19:
            case 20:
            case 21:
                {
                alt1=4;
                }
                break;
            case 22:
            case 23:
            case 24:
                {
                alt1=5;
                }
                break;
            case 25:
            case 26:
            case 27:
                {
                alt1=6;
                }
                break;
            case 28:
            case 29:
            case 30:
            case 31:
            case 32:
            case 33:
                {
                alt1=7;
                }
                break;
            case RULE_STRING:
            case RULE_ID:
                {
                alt1=8;
                }
                break;
            default:
                NoViableAltException nvae =
                    new NoViableAltException("", 1, 0, input);

                throw nvae;
            }

            switch (alt1) {
                case 1 :
                    // InternalElkGraphJson.g:967:2: ( ruleElkId )
                    {
                    // InternalElkGraphJson.g:967:2: ( ruleElkId )
                    // InternalElkGraphJson.g:968:3: ruleElkId
                    {
                     before(grammarAccess.getNodeElementAccess().getElkIdParserRuleCall_0()); 
                    pushFollow(FOLLOW_2);
                    ruleElkId();

                    state._fsp--;

                     after(grammarAccess.getNodeElementAccess().getElkIdParserRuleCall_0()); 

                    }


                    }
                    break;
                case 2 :
                    // InternalElkGraphJson.g:973:2: ( ruleShapeElement )
                    {
                    // InternalElkGraphJson.g:973:2: ( ruleShapeElement )
                    // InternalElkGraphJson.g:974:3: ruleShapeElement
                    {
                     before(grammarAccess.getNodeElementAccess().getShapeElementParserRuleCall_1()); 
                    pushFollow(FOLLOW_2);
                    ruleShapeElement();

                    state._fsp--;

                     after(grammarAccess.getNodeElementAccess().getShapeElementParserRuleCall_1()); 

                    }


                    }
                    break;
                case 3 :
                    // InternalElkGraphJson.g:979:2: ( ( rule__NodeElement__Group_2__0 ) )
                    {
                    // InternalElkGraphJson.g:979:2: ( ( rule__NodeElement__Group_2__0 ) )
                    // InternalElkGraphJson.g:980:3: ( rule__NodeElement__Group_2__0 )
                    {
                     before(grammarAccess.getNodeElementAccess().getGroup_2()); 
                    // InternalElkGraphJson.g:981:3: ( rule__NodeElement__Group_2__0 )
                    // InternalElkGraphJson.g:981:4: rule__NodeElement__Group_2__0
                    {
                    pushFollow(FOLLOW_2);
                    rule__NodeElement__Group_2__0();

                    state._fsp--;


                    }

                     after(grammarAccess.getNodeElementAccess().getGroup_2()); 

                    }


                    }
                    break;
                case 4 :
                    // InternalElkGraphJson.g:985:2: ( ( rule__NodeElement__Group_3__0 ) )
                    {
                    // InternalElkGraphJson.g:985:2: ( ( rule__NodeElement__Group_3__0 ) )
                    // InternalElkGraphJson.g:986:3: ( rule__NodeElement__Group_3__0 )
                    {
                     before(grammarAccess.getNodeElementAccess().getGroup_3()); 
                    // InternalElkGraphJson.g:987:3: ( rule__NodeElement__Group_3__0 )
                    // InternalElkGraphJson.g:987:4: rule__NodeElement__Group_3__0
                    {
                    pushFollow(FOLLOW_2);
                    rule__NodeElement__Group_3__0();

                    state._fsp--;


                    }

                     after(grammarAccess.getNodeElementAccess().getGroup_3()); 

                    }


                    }
                    break;
                case 5 :
                    // InternalElkGraphJson.g:991:2: ( ( rule__NodeElement__Group_4__0 ) )
                    {
                    // InternalElkGraphJson.g:991:2: ( ( rule__NodeElement__Group_4__0 ) )
                    // InternalElkGraphJson.g:992:3: ( rule__NodeElement__Group_4__0 )
                    {
                     before(grammarAccess.getNodeElementAccess().getGroup_4()); 
                    // InternalElkGraphJson.g:993:3: ( rule__NodeElement__Group_4__0 )
                    // InternalElkGraphJson.g:993:4: rule__NodeElement__Group_4__0
                    {
                    pushFollow(FOLLOW_2);
                    rule__NodeElement__Group_4__0();

                    state._fsp--;


                    }

                     after(grammarAccess.getNodeElementAccess().getGroup_4()); 

                    }


                    }
                    break;
                case 6 :
                    // InternalElkGraphJson.g:997:2: ( ( rule__NodeElement__Group_5__0 ) )
                    {
                    // InternalElkGraphJson.g:997:2: ( ( rule__NodeElement__Group_5__0 ) )
                    // InternalElkGraphJson.g:998:3: ( rule__NodeElement__Group_5__0 )
                    {
                     before(grammarAccess.getNodeElementAccess().getGroup_5()); 
                    // InternalElkGraphJson.g:999:3: ( rule__NodeElement__Group_5__0 )
                    // InternalElkGraphJson.g:999:4: rule__NodeElement__Group_5__0
                    {
                    pushFollow(FOLLOW_2);
                    rule__NodeElement__Group_5__0();

                    state._fsp--;


                    }

                     after(grammarAccess.getNodeElementAccess().getGroup_5()); 

                    }


                    }
                    break;
                case 7 :
                    // InternalElkGraphJson.g:1003:2: ( ( rule__NodeElement__Group_6__0 ) )
                    {
                    // InternalElkGraphJson.g:1003:2: ( ( rule__NodeElement__Group_6__0 ) )
                    // InternalElkGraphJson.g:1004:3: ( rule__NodeElement__Group_6__0 )
                    {
                     before(grammarAccess.getNodeElementAccess().getGroup_6()); 
                    // InternalElkGraphJson.g:1005:3: ( rule__NodeElement__Group_6__0 )
                    // InternalElkGraphJson.g:1005:4: rule__NodeElement__Group_6__0
                    {
                    pushFollow(FOLLOW_2);
                    rule__NodeElement__Group_6__0();

                    state._fsp--;


                    }

                     after(grammarAccess.getNodeElementAccess().getGroup_6()); 

                    }


                    }
                    break;
                case 8 :
                    // InternalElkGraphJson.g:1009:2: ( ruleJsonMember )
                    {
                    // InternalElkGraphJson.g:1009:2: ( ruleJsonMember )
                    // InternalElkGraphJson.g:1010:3: ruleJsonMember
                    {
                     before(grammarAccess.getNodeElementAccess().getJsonMemberParserRuleCall_7()); 
                    pushFollow(FOLLOW_2);
                    ruleJsonMember();

                    state._fsp--;

                     after(grammarAccess.getNodeElementAccess().getJsonMemberParserRuleCall_7()); 

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
    // $ANTLR end "rule__NodeElement__Alternatives"


    // $ANTLR start "rule__PortElement__Alternatives"
    // InternalElkGraphJson.g:1019:1: rule__PortElement__Alternatives : ( ( ruleElkId ) | ( ruleShapeElement ) | ( ( rule__PortElement__Group_2__0 ) ) | ( ( rule__PortElement__Group_3__0 ) ) | ( ruleJsonMember ) );
    public final void rule__PortElement__Alternatives() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:1023:1: ( ( ruleElkId ) | ( ruleShapeElement ) | ( ( rule__PortElement__Group_2__0 ) ) | ( ( rule__PortElement__Group_3__0 ) ) | ( ruleJsonMember ) )
            int alt2=5;
            switch ( input.LA(1) ) {
            case 34:
            case 35:
            case 36:
                {
                alt2=1;
                }
                break;
            case 37:
            case 38:
            case 39:
            case 40:
            case 41:
            case 42:
            case 43:
            case 44:
            case 45:
            case 46:
            case 47:
            case 48:
                {
                alt2=2;
                }
                break;
            case 22:
            case 23:
            case 24:
                {
                alt2=3;
                }
                break;
            case 28:
            case 29:
            case 30:
            case 31:
            case 32:
            case 33:
                {
                alt2=4;
                }
                break;
            case RULE_STRING:
            case RULE_ID:
                {
                alt2=5;
                }
                break;
            default:
                NoViableAltException nvae =
                    new NoViableAltException("", 2, 0, input);

                throw nvae;
            }

            switch (alt2) {
                case 1 :
                    // InternalElkGraphJson.g:1024:2: ( ruleElkId )
                    {
                    // InternalElkGraphJson.g:1024:2: ( ruleElkId )
                    // InternalElkGraphJson.g:1025:3: ruleElkId
                    {
                     before(grammarAccess.getPortElementAccess().getElkIdParserRuleCall_0()); 
                    pushFollow(FOLLOW_2);
                    ruleElkId();

                    state._fsp--;

                     after(grammarAccess.getPortElementAccess().getElkIdParserRuleCall_0()); 

                    }


                    }
                    break;
                case 2 :
                    // InternalElkGraphJson.g:1030:2: ( ruleShapeElement )
                    {
                    // InternalElkGraphJson.g:1030:2: ( ruleShapeElement )
                    // InternalElkGraphJson.g:1031:3: ruleShapeElement
                    {
                     before(grammarAccess.getPortElementAccess().getShapeElementParserRuleCall_1()); 
                    pushFollow(FOLLOW_2);
                    ruleShapeElement();

                    state._fsp--;

                     after(grammarAccess.getPortElementAccess().getShapeElementParserRuleCall_1()); 

                    }


                    }
                    break;
                case 3 :
                    // InternalElkGraphJson.g:1036:2: ( ( rule__PortElement__Group_2__0 ) )
                    {
                    // InternalElkGraphJson.g:1036:2: ( ( rule__PortElement__Group_2__0 ) )
                    // InternalElkGraphJson.g:1037:3: ( rule__PortElement__Group_2__0 )
                    {
                     before(grammarAccess.getPortElementAccess().getGroup_2()); 
                    // InternalElkGraphJson.g:1038:3: ( rule__PortElement__Group_2__0 )
                    // InternalElkGraphJson.g:1038:4: rule__PortElement__Group_2__0
                    {
                    pushFollow(FOLLOW_2);
                    rule__PortElement__Group_2__0();

                    state._fsp--;


                    }

                     after(grammarAccess.getPortElementAccess().getGroup_2()); 

                    }


                    }
                    break;
                case 4 :
                    // InternalElkGraphJson.g:1042:2: ( ( rule__PortElement__Group_3__0 ) )
                    {
                    // InternalElkGraphJson.g:1042:2: ( ( rule__PortElement__Group_3__0 ) )
                    // InternalElkGraphJson.g:1043:3: ( rule__PortElement__Group_3__0 )
                    {
                     before(grammarAccess.getPortElementAccess().getGroup_3()); 
                    // InternalElkGraphJson.g:1044:3: ( rule__PortElement__Group_3__0 )
                    // InternalElkGraphJson.g:1044:4: rule__PortElement__Group_3__0
                    {
                    pushFollow(FOLLOW_2);
                    rule__PortElement__Group_3__0();

                    state._fsp--;


                    }

                     after(grammarAccess.getPortElementAccess().getGroup_3()); 

                    }


                    }
                    break;
                case 5 :
                    // InternalElkGraphJson.g:1048:2: ( ruleJsonMember )
                    {
                    // InternalElkGraphJson.g:1048:2: ( ruleJsonMember )
                    // InternalElkGraphJson.g:1049:3: ruleJsonMember
                    {
                     before(grammarAccess.getPortElementAccess().getJsonMemberParserRuleCall_4()); 
                    pushFollow(FOLLOW_2);
                    ruleJsonMember();

                    state._fsp--;

                     after(grammarAccess.getPortElementAccess().getJsonMemberParserRuleCall_4()); 

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
    // $ANTLR end "rule__PortElement__Alternatives"


    // $ANTLR start "rule__LabelElement__Alternatives"
    // InternalElkGraphJson.g:1058:1: rule__LabelElement__Alternatives : ( ( ruleElkId ) | ( ruleShapeElement ) | ( ( rule__LabelElement__Group_2__0 ) ) | ( ( rule__LabelElement__Group_3__0 ) ) | ( ( rule__LabelElement__Group_4__0 ) ) | ( ruleJsonMember ) );
    public final void rule__LabelElement__Alternatives() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:1062:1: ( ( ruleElkId ) | ( ruleShapeElement ) | ( ( rule__LabelElement__Group_2__0 ) ) | ( ( rule__LabelElement__Group_3__0 ) ) | ( ( rule__LabelElement__Group_4__0 ) ) | ( ruleJsonMember ) )
            int alt3=6;
            switch ( input.LA(1) ) {
            case 34:
            case 35:
            case 36:
                {
                alt3=1;
                }
                break;
            case 37:
            case 38:
            case 39:
            case 40:
            case 41:
            case 42:
            case 43:
            case 44:
            case 45:
            case 46:
            case 47:
            case 48:
                {
                alt3=2;
                }
                break;
            case 55:
            case 56:
            case 57:
                {
                alt3=3;
                }
                break;
            case 22:
            case 23:
            case 24:
                {
                alt3=4;
                }
                break;
            case 28:
            case 29:
            case 30:
            case 31:
            case 32:
            case 33:
                {
                alt3=5;
                }
                break;
            case RULE_STRING:
            case RULE_ID:
                {
                alt3=6;
                }
                break;
            default:
                NoViableAltException nvae =
                    new NoViableAltException("", 3, 0, input);

                throw nvae;
            }

            switch (alt3) {
                case 1 :
                    // InternalElkGraphJson.g:1063:2: ( ruleElkId )
                    {
                    // InternalElkGraphJson.g:1063:2: ( ruleElkId )
                    // InternalElkGraphJson.g:1064:3: ruleElkId
                    {
                     before(grammarAccess.getLabelElementAccess().getElkIdParserRuleCall_0()); 
                    pushFollow(FOLLOW_2);
                    ruleElkId();

                    state._fsp--;

                     after(grammarAccess.getLabelElementAccess().getElkIdParserRuleCall_0()); 

                    }


                    }
                    break;
                case 2 :
                    // InternalElkGraphJson.g:1069:2: ( ruleShapeElement )
                    {
                    // InternalElkGraphJson.g:1069:2: ( ruleShapeElement )
                    // InternalElkGraphJson.g:1070:3: ruleShapeElement
                    {
                     before(grammarAccess.getLabelElementAccess().getShapeElementParserRuleCall_1()); 
                    pushFollow(FOLLOW_2);
                    ruleShapeElement();

                    state._fsp--;

                     after(grammarAccess.getLabelElementAccess().getShapeElementParserRuleCall_1()); 

                    }


                    }
                    break;
                case 3 :
                    // InternalElkGraphJson.g:1075:2: ( ( rule__LabelElement__Group_2__0 ) )
                    {
                    // InternalElkGraphJson.g:1075:2: ( ( rule__LabelElement__Group_2__0 ) )
                    // InternalElkGraphJson.g:1076:3: ( rule__LabelElement__Group_2__0 )
                    {
                     before(grammarAccess.getLabelElementAccess().getGroup_2()); 
                    // InternalElkGraphJson.g:1077:3: ( rule__LabelElement__Group_2__0 )
                    // InternalElkGraphJson.g:1077:4: rule__LabelElement__Group_2__0
                    {
                    pushFollow(FOLLOW_2);
                    rule__LabelElement__Group_2__0();

                    state._fsp--;


                    }

                     after(grammarAccess.getLabelElementAccess().getGroup_2()); 

                    }


                    }
                    break;
                case 4 :
                    // InternalElkGraphJson.g:1081:2: ( ( rule__LabelElement__Group_3__0 ) )
                    {
                    // InternalElkGraphJson.g:1081:2: ( ( rule__LabelElement__Group_3__0 ) )
                    // InternalElkGraphJson.g:1082:3: ( rule__LabelElement__Group_3__0 )
                    {
                     before(grammarAccess.getLabelElementAccess().getGroup_3()); 
                    // InternalElkGraphJson.g:1083:3: ( rule__LabelElement__Group_3__0 )
                    // InternalElkGraphJson.g:1083:4: rule__LabelElement__Group_3__0
                    {
                    pushFollow(FOLLOW_2);
                    rule__LabelElement__Group_3__0();

                    state._fsp--;


                    }

                     after(grammarAccess.getLabelElementAccess().getGroup_3()); 

                    }


                    }
                    break;
                case 5 :
                    // InternalElkGraphJson.g:1087:2: ( ( rule__LabelElement__Group_4__0 ) )
                    {
                    // InternalElkGraphJson.g:1087:2: ( ( rule__LabelElement__Group_4__0 ) )
                    // InternalElkGraphJson.g:1088:3: ( rule__LabelElement__Group_4__0 )
                    {
                     before(grammarAccess.getLabelElementAccess().getGroup_4()); 
                    // InternalElkGraphJson.g:1089:3: ( rule__LabelElement__Group_4__0 )
                    // InternalElkGraphJson.g:1089:4: rule__LabelElement__Group_4__0
                    {
                    pushFollow(FOLLOW_2);
                    rule__LabelElement__Group_4__0();

                    state._fsp--;


                    }

                     after(grammarAccess.getLabelElementAccess().getGroup_4()); 

                    }


                    }
                    break;
                case 6 :
                    // InternalElkGraphJson.g:1093:2: ( ruleJsonMember )
                    {
                    // InternalElkGraphJson.g:1093:2: ( ruleJsonMember )
                    // InternalElkGraphJson.g:1094:3: ruleJsonMember
                    {
                     before(grammarAccess.getLabelElementAccess().getJsonMemberParserRuleCall_5()); 
                    pushFollow(FOLLOW_2);
                    ruleJsonMember();

                    state._fsp--;

                     after(grammarAccess.getLabelElementAccess().getJsonMemberParserRuleCall_5()); 

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
    // $ANTLR end "rule__LabelElement__Alternatives"


    // $ANTLR start "rule__EdgeElement__Alternatives"
    // InternalElkGraphJson.g:1103:1: rule__EdgeElement__Alternatives : ( ( ruleElkId ) | ( ( rule__EdgeElement__Group_1__0 ) ) | ( ( rule__EdgeElement__Group_2__0 ) ) | ( ( rule__EdgeElement__Group_3__0 ) ) | ( ( rule__EdgeElement__Group_4__0 ) ) | ( ruleJsonMember ) );
    public final void rule__EdgeElement__Alternatives() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:1107:1: ( ( ruleElkId ) | ( ( rule__EdgeElement__Group_1__0 ) ) | ( ( rule__EdgeElement__Group_2__0 ) ) | ( ( rule__EdgeElement__Group_3__0 ) ) | ( ( rule__EdgeElement__Group_4__0 ) ) | ( ruleJsonMember ) )
            int alt4=6;
            switch ( input.LA(1) ) {
            case 34:
            case 35:
            case 36:
                {
                alt4=1;
                }
                break;
            case 49:
            case 50:
            case 51:
                {
                alt4=2;
                }
                break;
            case 52:
            case 53:
            case 54:
                {
                alt4=3;
                }
                break;
            case 22:
            case 23:
            case 24:
                {
                alt4=4;
                }
                break;
            case 28:
            case 29:
            case 30:
            case 31:
            case 32:
            case 33:
                {
                alt4=5;
                }
                break;
            case RULE_STRING:
            case RULE_ID:
                {
                alt4=6;
                }
                break;
            default:
                NoViableAltException nvae =
                    new NoViableAltException("", 4, 0, input);

                throw nvae;
            }

            switch (alt4) {
                case 1 :
                    // InternalElkGraphJson.g:1108:2: ( ruleElkId )
                    {
                    // InternalElkGraphJson.g:1108:2: ( ruleElkId )
                    // InternalElkGraphJson.g:1109:3: ruleElkId
                    {
                     before(grammarAccess.getEdgeElementAccess().getElkIdParserRuleCall_0()); 
                    pushFollow(FOLLOW_2);
                    ruleElkId();

                    state._fsp--;

                     after(grammarAccess.getEdgeElementAccess().getElkIdParserRuleCall_0()); 

                    }


                    }
                    break;
                case 2 :
                    // InternalElkGraphJson.g:1114:2: ( ( rule__EdgeElement__Group_1__0 ) )
                    {
                    // InternalElkGraphJson.g:1114:2: ( ( rule__EdgeElement__Group_1__0 ) )
                    // InternalElkGraphJson.g:1115:3: ( rule__EdgeElement__Group_1__0 )
                    {
                     before(grammarAccess.getEdgeElementAccess().getGroup_1()); 
                    // InternalElkGraphJson.g:1116:3: ( rule__EdgeElement__Group_1__0 )
                    // InternalElkGraphJson.g:1116:4: rule__EdgeElement__Group_1__0
                    {
                    pushFollow(FOLLOW_2);
                    rule__EdgeElement__Group_1__0();

                    state._fsp--;


                    }

                     after(grammarAccess.getEdgeElementAccess().getGroup_1()); 

                    }


                    }
                    break;
                case 3 :
                    // InternalElkGraphJson.g:1120:2: ( ( rule__EdgeElement__Group_2__0 ) )
                    {
                    // InternalElkGraphJson.g:1120:2: ( ( rule__EdgeElement__Group_2__0 ) )
                    // InternalElkGraphJson.g:1121:3: ( rule__EdgeElement__Group_2__0 )
                    {
                     before(grammarAccess.getEdgeElementAccess().getGroup_2()); 
                    // InternalElkGraphJson.g:1122:3: ( rule__EdgeElement__Group_2__0 )
                    // InternalElkGraphJson.g:1122:4: rule__EdgeElement__Group_2__0
                    {
                    pushFollow(FOLLOW_2);
                    rule__EdgeElement__Group_2__0();

                    state._fsp--;


                    }

                     after(grammarAccess.getEdgeElementAccess().getGroup_2()); 

                    }


                    }
                    break;
                case 4 :
                    // InternalElkGraphJson.g:1126:2: ( ( rule__EdgeElement__Group_3__0 ) )
                    {
                    // InternalElkGraphJson.g:1126:2: ( ( rule__EdgeElement__Group_3__0 ) )
                    // InternalElkGraphJson.g:1127:3: ( rule__EdgeElement__Group_3__0 )
                    {
                     before(grammarAccess.getEdgeElementAccess().getGroup_3()); 
                    // InternalElkGraphJson.g:1128:3: ( rule__EdgeElement__Group_3__0 )
                    // InternalElkGraphJson.g:1128:4: rule__EdgeElement__Group_3__0
                    {
                    pushFollow(FOLLOW_2);
                    rule__EdgeElement__Group_3__0();

                    state._fsp--;


                    }

                     after(grammarAccess.getEdgeElementAccess().getGroup_3()); 

                    }


                    }
                    break;
                case 5 :
                    // InternalElkGraphJson.g:1132:2: ( ( rule__EdgeElement__Group_4__0 ) )
                    {
                    // InternalElkGraphJson.g:1132:2: ( ( rule__EdgeElement__Group_4__0 ) )
                    // InternalElkGraphJson.g:1133:3: ( rule__EdgeElement__Group_4__0 )
                    {
                     before(grammarAccess.getEdgeElementAccess().getGroup_4()); 
                    // InternalElkGraphJson.g:1134:3: ( rule__EdgeElement__Group_4__0 )
                    // InternalElkGraphJson.g:1134:4: rule__EdgeElement__Group_4__0
                    {
                    pushFollow(FOLLOW_2);
                    rule__EdgeElement__Group_4__0();

                    state._fsp--;


                    }

                     after(grammarAccess.getEdgeElementAccess().getGroup_4()); 

                    }


                    }
                    break;
                case 6 :
                    // InternalElkGraphJson.g:1138:2: ( ruleJsonMember )
                    {
                    // InternalElkGraphJson.g:1138:2: ( ruleJsonMember )
                    // InternalElkGraphJson.g:1139:3: ruleJsonMember
                    {
                     before(grammarAccess.getEdgeElementAccess().getJsonMemberParserRuleCall_5()); 
                    pushFollow(FOLLOW_2);
                    ruleJsonMember();

                    state._fsp--;

                     after(grammarAccess.getEdgeElementAccess().getJsonMemberParserRuleCall_5()); 

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
    // $ANTLR end "rule__EdgeElement__Alternatives"


    // $ANTLR start "rule__ShapeElement__Alternatives"
    // InternalElkGraphJson.g:1148:1: rule__ShapeElement__Alternatives : ( ( ( rule__ShapeElement__Group_0__0 ) ) | ( ( rule__ShapeElement__Group_1__0 ) ) | ( ( rule__ShapeElement__Group_2__0 ) ) | ( ( rule__ShapeElement__Group_3__0 ) ) );
    public final void rule__ShapeElement__Alternatives() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:1152:1: ( ( ( rule__ShapeElement__Group_0__0 ) ) | ( ( rule__ShapeElement__Group_1__0 ) ) | ( ( rule__ShapeElement__Group_2__0 ) ) | ( ( rule__ShapeElement__Group_3__0 ) ) )
            int alt5=4;
            switch ( input.LA(1) ) {
            case 37:
            case 38:
            case 39:
                {
                alt5=1;
                }
                break;
            case 40:
            case 41:
            case 42:
                {
                alt5=2;
                }
                break;
            case 43:
            case 44:
            case 45:
                {
                alt5=3;
                }
                break;
            case 46:
            case 47:
            case 48:
                {
                alt5=4;
                }
                break;
            default:
                NoViableAltException nvae =
                    new NoViableAltException("", 5, 0, input);

                throw nvae;
            }

            switch (alt5) {
                case 1 :
                    // InternalElkGraphJson.g:1153:2: ( ( rule__ShapeElement__Group_0__0 ) )
                    {
                    // InternalElkGraphJson.g:1153:2: ( ( rule__ShapeElement__Group_0__0 ) )
                    // InternalElkGraphJson.g:1154:3: ( rule__ShapeElement__Group_0__0 )
                    {
                     before(grammarAccess.getShapeElementAccess().getGroup_0()); 
                    // InternalElkGraphJson.g:1155:3: ( rule__ShapeElement__Group_0__0 )
                    // InternalElkGraphJson.g:1155:4: rule__ShapeElement__Group_0__0
                    {
                    pushFollow(FOLLOW_2);
                    rule__ShapeElement__Group_0__0();

                    state._fsp--;


                    }

                     after(grammarAccess.getShapeElementAccess().getGroup_0()); 

                    }


                    }
                    break;
                case 2 :
                    // InternalElkGraphJson.g:1159:2: ( ( rule__ShapeElement__Group_1__0 ) )
                    {
                    // InternalElkGraphJson.g:1159:2: ( ( rule__ShapeElement__Group_1__0 ) )
                    // InternalElkGraphJson.g:1160:3: ( rule__ShapeElement__Group_1__0 )
                    {
                     before(grammarAccess.getShapeElementAccess().getGroup_1()); 
                    // InternalElkGraphJson.g:1161:3: ( rule__ShapeElement__Group_1__0 )
                    // InternalElkGraphJson.g:1161:4: rule__ShapeElement__Group_1__0
                    {
                    pushFollow(FOLLOW_2);
                    rule__ShapeElement__Group_1__0();

                    state._fsp--;


                    }

                     after(grammarAccess.getShapeElementAccess().getGroup_1()); 

                    }


                    }
                    break;
                case 3 :
                    // InternalElkGraphJson.g:1165:2: ( ( rule__ShapeElement__Group_2__0 ) )
                    {
                    // InternalElkGraphJson.g:1165:2: ( ( rule__ShapeElement__Group_2__0 ) )
                    // InternalElkGraphJson.g:1166:3: ( rule__ShapeElement__Group_2__0 )
                    {
                     before(grammarAccess.getShapeElementAccess().getGroup_2()); 
                    // InternalElkGraphJson.g:1167:3: ( rule__ShapeElement__Group_2__0 )
                    // InternalElkGraphJson.g:1167:4: rule__ShapeElement__Group_2__0
                    {
                    pushFollow(FOLLOW_2);
                    rule__ShapeElement__Group_2__0();

                    state._fsp--;


                    }

                     after(grammarAccess.getShapeElementAccess().getGroup_2()); 

                    }


                    }
                    break;
                case 4 :
                    // InternalElkGraphJson.g:1171:2: ( ( rule__ShapeElement__Group_3__0 ) )
                    {
                    // InternalElkGraphJson.g:1171:2: ( ( rule__ShapeElement__Group_3__0 ) )
                    // InternalElkGraphJson.g:1172:3: ( rule__ShapeElement__Group_3__0 )
                    {
                     before(grammarAccess.getShapeElementAccess().getGroup_3()); 
                    // InternalElkGraphJson.g:1173:3: ( rule__ShapeElement__Group_3__0 )
                    // InternalElkGraphJson.g:1173:4: rule__ShapeElement__Group_3__0
                    {
                    pushFollow(FOLLOW_2);
                    rule__ShapeElement__Group_3__0();

                    state._fsp--;


                    }

                     after(grammarAccess.getShapeElementAccess().getGroup_3()); 

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
    // $ANTLR end "rule__ShapeElement__Alternatives"


    // $ANTLR start "rule__Property__Alternatives_2"
    // InternalElkGraphJson.g:1181:1: rule__Property__Alternatives_2 : ( ( ( rule__Property__ValueAssignment_2_0 ) ) | ( ( rule__Property__ValueAssignment_2_1 ) ) | ( ( rule__Property__ValueAssignment_2_2 ) ) );
    public final void rule__Property__Alternatives_2() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:1185:1: ( ( ( rule__Property__ValueAssignment_2_0 ) ) | ( ( rule__Property__ValueAssignment_2_1 ) ) | ( ( rule__Property__ValueAssignment_2_2 ) ) )
            int alt6=3;
            switch ( input.LA(1) ) {
            case RULE_STRING:
                {
                alt6=1;
                }
                break;
            case RULE_SIGNED_INT:
            case RULE_FLOAT:
                {
                alt6=2;
                }
                break;
            case 13:
            case 14:
                {
                alt6=3;
                }
                break;
            default:
                NoViableAltException nvae =
                    new NoViableAltException("", 6, 0, input);

                throw nvae;
            }

            switch (alt6) {
                case 1 :
                    // InternalElkGraphJson.g:1186:2: ( ( rule__Property__ValueAssignment_2_0 ) )
                    {
                    // InternalElkGraphJson.g:1186:2: ( ( rule__Property__ValueAssignment_2_0 ) )
                    // InternalElkGraphJson.g:1187:3: ( rule__Property__ValueAssignment_2_0 )
                    {
                     before(grammarAccess.getPropertyAccess().getValueAssignment_2_0()); 
                    // InternalElkGraphJson.g:1188:3: ( rule__Property__ValueAssignment_2_0 )
                    // InternalElkGraphJson.g:1188:4: rule__Property__ValueAssignment_2_0
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
                    // InternalElkGraphJson.g:1192:2: ( ( rule__Property__ValueAssignment_2_1 ) )
                    {
                    // InternalElkGraphJson.g:1192:2: ( ( rule__Property__ValueAssignment_2_1 ) )
                    // InternalElkGraphJson.g:1193:3: ( rule__Property__ValueAssignment_2_1 )
                    {
                     before(grammarAccess.getPropertyAccess().getValueAssignment_2_1()); 
                    // InternalElkGraphJson.g:1194:3: ( rule__Property__ValueAssignment_2_1 )
                    // InternalElkGraphJson.g:1194:4: rule__Property__ValueAssignment_2_1
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
                    // InternalElkGraphJson.g:1198:2: ( ( rule__Property__ValueAssignment_2_2 ) )
                    {
                    // InternalElkGraphJson.g:1198:2: ( ( rule__Property__ValueAssignment_2_2 ) )
                    // InternalElkGraphJson.g:1199:3: ( rule__Property__ValueAssignment_2_2 )
                    {
                     before(grammarAccess.getPropertyAccess().getValueAssignment_2_2()); 
                    // InternalElkGraphJson.g:1200:3: ( rule__Property__ValueAssignment_2_2 )
                    // InternalElkGraphJson.g:1200:4: rule__Property__ValueAssignment_2_2
                    {
                    pushFollow(FOLLOW_2);
                    rule__Property__ValueAssignment_2_2();

                    state._fsp--;


                    }

                     after(grammarAccess.getPropertyAccess().getValueAssignment_2_2()); 

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


    // $ANTLR start "rule__PropertyKey__Alternatives"
    // InternalElkGraphJson.g:1208:1: rule__PropertyKey__Alternatives : ( ( RULE_STRING ) | ( RULE_ID ) );
    public final void rule__PropertyKey__Alternatives() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:1212:1: ( ( RULE_STRING ) | ( RULE_ID ) )
            int alt7=2;
            int LA7_0 = input.LA(1);

            if ( (LA7_0==RULE_STRING) ) {
                alt7=1;
            }
            else if ( (LA7_0==RULE_ID) ) {
                alt7=2;
            }
            else {
                NoViableAltException nvae =
                    new NoViableAltException("", 7, 0, input);

                throw nvae;
            }
            switch (alt7) {
                case 1 :
                    // InternalElkGraphJson.g:1213:2: ( RULE_STRING )
                    {
                    // InternalElkGraphJson.g:1213:2: ( RULE_STRING )
                    // InternalElkGraphJson.g:1214:3: RULE_STRING
                    {
                     before(grammarAccess.getPropertyKeyAccess().getSTRINGTerminalRuleCall_0()); 
                    match(input,RULE_STRING,FOLLOW_2); 
                     after(grammarAccess.getPropertyKeyAccess().getSTRINGTerminalRuleCall_0()); 

                    }


                    }
                    break;
                case 2 :
                    // InternalElkGraphJson.g:1219:2: ( RULE_ID )
                    {
                    // InternalElkGraphJson.g:1219:2: ( RULE_ID )
                    // InternalElkGraphJson.g:1220:3: RULE_ID
                    {
                     before(grammarAccess.getPropertyKeyAccess().getIDTerminalRuleCall_1()); 
                    match(input,RULE_ID,FOLLOW_2); 
                     after(grammarAccess.getPropertyKeyAccess().getIDTerminalRuleCall_1()); 

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
    // $ANTLR end "rule__PropertyKey__Alternatives"


    // $ANTLR start "rule__NumberValue__Alternatives"
    // InternalElkGraphJson.g:1229:1: rule__NumberValue__Alternatives : ( ( RULE_SIGNED_INT ) | ( RULE_FLOAT ) );
    public final void rule__NumberValue__Alternatives() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:1233:1: ( ( RULE_SIGNED_INT ) | ( RULE_FLOAT ) )
            int alt8=2;
            int LA8_0 = input.LA(1);

            if ( (LA8_0==RULE_SIGNED_INT) ) {
                alt8=1;
            }
            else if ( (LA8_0==RULE_FLOAT) ) {
                alt8=2;
            }
            else {
                NoViableAltException nvae =
                    new NoViableAltException("", 8, 0, input);

                throw nvae;
            }
            switch (alt8) {
                case 1 :
                    // InternalElkGraphJson.g:1234:2: ( RULE_SIGNED_INT )
                    {
                    // InternalElkGraphJson.g:1234:2: ( RULE_SIGNED_INT )
                    // InternalElkGraphJson.g:1235:3: RULE_SIGNED_INT
                    {
                     before(grammarAccess.getNumberValueAccess().getSIGNED_INTTerminalRuleCall_0()); 
                    match(input,RULE_SIGNED_INT,FOLLOW_2); 
                     after(grammarAccess.getNumberValueAccess().getSIGNED_INTTerminalRuleCall_0()); 

                    }


                    }
                    break;
                case 2 :
                    // InternalElkGraphJson.g:1240:2: ( RULE_FLOAT )
                    {
                    // InternalElkGraphJson.g:1240:2: ( RULE_FLOAT )
                    // InternalElkGraphJson.g:1241:3: RULE_FLOAT
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
    // InternalElkGraphJson.g:1250:1: rule__BooleanValue__Alternatives : ( ( 'true' ) | ( 'false' ) );
    public final void rule__BooleanValue__Alternatives() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:1254:1: ( ( 'true' ) | ( 'false' ) )
            int alt9=2;
            int LA9_0 = input.LA(1);

            if ( (LA9_0==13) ) {
                alt9=1;
            }
            else if ( (LA9_0==14) ) {
                alt9=2;
            }
            else {
                NoViableAltException nvae =
                    new NoViableAltException("", 9, 0, input);

                throw nvae;
            }
            switch (alt9) {
                case 1 :
                    // InternalElkGraphJson.g:1255:2: ( 'true' )
                    {
                    // InternalElkGraphJson.g:1255:2: ( 'true' )
                    // InternalElkGraphJson.g:1256:3: 'true'
                    {
                     before(grammarAccess.getBooleanValueAccess().getTrueKeyword_0()); 
                    match(input,13,FOLLOW_2); 
                     after(grammarAccess.getBooleanValueAccess().getTrueKeyword_0()); 

                    }


                    }
                    break;
                case 2 :
                    // InternalElkGraphJson.g:1261:2: ( 'false' )
                    {
                    // InternalElkGraphJson.g:1261:2: ( 'false' )
                    // InternalElkGraphJson.g:1262:3: 'false'
                    {
                     before(grammarAccess.getBooleanValueAccess().getFalseKeyword_1()); 
                    match(input,14,FOLLOW_2); 
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


    // $ANTLR start "rule__Number__Alternatives"
    // InternalElkGraphJson.g:1271:1: rule__Number__Alternatives : ( ( RULE_SIGNED_INT ) | ( RULE_FLOAT ) );
    public final void rule__Number__Alternatives() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:1275:1: ( ( RULE_SIGNED_INT ) | ( RULE_FLOAT ) )
            int alt10=2;
            int LA10_0 = input.LA(1);

            if ( (LA10_0==RULE_SIGNED_INT) ) {
                alt10=1;
            }
            else if ( (LA10_0==RULE_FLOAT) ) {
                alt10=2;
            }
            else {
                NoViableAltException nvae =
                    new NoViableAltException("", 10, 0, input);

                throw nvae;
            }
            switch (alt10) {
                case 1 :
                    // InternalElkGraphJson.g:1276:2: ( RULE_SIGNED_INT )
                    {
                    // InternalElkGraphJson.g:1276:2: ( RULE_SIGNED_INT )
                    // InternalElkGraphJson.g:1277:3: RULE_SIGNED_INT
                    {
                     before(grammarAccess.getNumberAccess().getSIGNED_INTTerminalRuleCall_0()); 
                    match(input,RULE_SIGNED_INT,FOLLOW_2); 
                     after(grammarAccess.getNumberAccess().getSIGNED_INTTerminalRuleCall_0()); 

                    }


                    }
                    break;
                case 2 :
                    // InternalElkGraphJson.g:1282:2: ( RULE_FLOAT )
                    {
                    // InternalElkGraphJson.g:1282:2: ( RULE_FLOAT )
                    // InternalElkGraphJson.g:1283:3: RULE_FLOAT
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


    // $ANTLR start "rule__JsonMember__Alternatives_0"
    // InternalElkGraphJson.g:1292:1: rule__JsonMember__Alternatives_0 : ( ( RULE_STRING ) | ( RULE_ID ) );
    public final void rule__JsonMember__Alternatives_0() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:1296:1: ( ( RULE_STRING ) | ( RULE_ID ) )
            int alt11=2;
            int LA11_0 = input.LA(1);

            if ( (LA11_0==RULE_STRING) ) {
                alt11=1;
            }
            else if ( (LA11_0==RULE_ID) ) {
                alt11=2;
            }
            else {
                NoViableAltException nvae =
                    new NoViableAltException("", 11, 0, input);

                throw nvae;
            }
            switch (alt11) {
                case 1 :
                    // InternalElkGraphJson.g:1297:2: ( RULE_STRING )
                    {
                    // InternalElkGraphJson.g:1297:2: ( RULE_STRING )
                    // InternalElkGraphJson.g:1298:3: RULE_STRING
                    {
                     before(grammarAccess.getJsonMemberAccess().getSTRINGTerminalRuleCall_0_0()); 
                    match(input,RULE_STRING,FOLLOW_2); 
                     after(grammarAccess.getJsonMemberAccess().getSTRINGTerminalRuleCall_0_0()); 

                    }


                    }
                    break;
                case 2 :
                    // InternalElkGraphJson.g:1303:2: ( RULE_ID )
                    {
                    // InternalElkGraphJson.g:1303:2: ( RULE_ID )
                    // InternalElkGraphJson.g:1304:3: RULE_ID
                    {
                     before(grammarAccess.getJsonMemberAccess().getIDTerminalRuleCall_0_1()); 
                    match(input,RULE_ID,FOLLOW_2); 
                     after(grammarAccess.getJsonMemberAccess().getIDTerminalRuleCall_0_1()); 

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
    // $ANTLR end "rule__JsonMember__Alternatives_0"


    // $ANTLR start "rule__JsonValue__Alternatives"
    // InternalElkGraphJson.g:1313:1: rule__JsonValue__Alternatives : ( ( RULE_STRING ) | ( ruleNumber ) | ( ruleJsonObject ) | ( ruleJsonArray ) | ( ruleBooleanValue ) | ( 'null' ) );
    public final void rule__JsonValue__Alternatives() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:1317:1: ( ( RULE_STRING ) | ( ruleNumber ) | ( ruleJsonObject ) | ( ruleJsonArray ) | ( ruleBooleanValue ) | ( 'null' ) )
            int alt12=6;
            switch ( input.LA(1) ) {
            case RULE_STRING:
                {
                alt12=1;
                }
                break;
            case RULE_SIGNED_INT:
            case RULE_FLOAT:
                {
                alt12=2;
                }
                break;
            case 58:
                {
                alt12=3;
                }
                break;
            case 62:
                {
                alt12=4;
                }
                break;
            case 13:
            case 14:
                {
                alt12=5;
                }
                break;
            case 15:
                {
                alt12=6;
                }
                break;
            default:
                NoViableAltException nvae =
                    new NoViableAltException("", 12, 0, input);

                throw nvae;
            }

            switch (alt12) {
                case 1 :
                    // InternalElkGraphJson.g:1318:2: ( RULE_STRING )
                    {
                    // InternalElkGraphJson.g:1318:2: ( RULE_STRING )
                    // InternalElkGraphJson.g:1319:3: RULE_STRING
                    {
                     before(grammarAccess.getJsonValueAccess().getSTRINGTerminalRuleCall_0()); 
                    match(input,RULE_STRING,FOLLOW_2); 
                     after(grammarAccess.getJsonValueAccess().getSTRINGTerminalRuleCall_0()); 

                    }


                    }
                    break;
                case 2 :
                    // InternalElkGraphJson.g:1324:2: ( ruleNumber )
                    {
                    // InternalElkGraphJson.g:1324:2: ( ruleNumber )
                    // InternalElkGraphJson.g:1325:3: ruleNumber
                    {
                     before(grammarAccess.getJsonValueAccess().getNumberParserRuleCall_1()); 
                    pushFollow(FOLLOW_2);
                    ruleNumber();

                    state._fsp--;

                     after(grammarAccess.getJsonValueAccess().getNumberParserRuleCall_1()); 

                    }


                    }
                    break;
                case 3 :
                    // InternalElkGraphJson.g:1330:2: ( ruleJsonObject )
                    {
                    // InternalElkGraphJson.g:1330:2: ( ruleJsonObject )
                    // InternalElkGraphJson.g:1331:3: ruleJsonObject
                    {
                     before(grammarAccess.getJsonValueAccess().getJsonObjectParserRuleCall_2()); 
                    pushFollow(FOLLOW_2);
                    ruleJsonObject();

                    state._fsp--;

                     after(grammarAccess.getJsonValueAccess().getJsonObjectParserRuleCall_2()); 

                    }


                    }
                    break;
                case 4 :
                    // InternalElkGraphJson.g:1336:2: ( ruleJsonArray )
                    {
                    // InternalElkGraphJson.g:1336:2: ( ruleJsonArray )
                    // InternalElkGraphJson.g:1337:3: ruleJsonArray
                    {
                     before(grammarAccess.getJsonValueAccess().getJsonArrayParserRuleCall_3()); 
                    pushFollow(FOLLOW_2);
                    ruleJsonArray();

                    state._fsp--;

                     after(grammarAccess.getJsonValueAccess().getJsonArrayParserRuleCall_3()); 

                    }


                    }
                    break;
                case 5 :
                    // InternalElkGraphJson.g:1342:2: ( ruleBooleanValue )
                    {
                    // InternalElkGraphJson.g:1342:2: ( ruleBooleanValue )
                    // InternalElkGraphJson.g:1343:3: ruleBooleanValue
                    {
                     before(grammarAccess.getJsonValueAccess().getBooleanValueParserRuleCall_4()); 
                    pushFollow(FOLLOW_2);
                    ruleBooleanValue();

                    state._fsp--;

                     after(grammarAccess.getJsonValueAccess().getBooleanValueParserRuleCall_4()); 

                    }


                    }
                    break;
                case 6 :
                    // InternalElkGraphJson.g:1348:2: ( 'null' )
                    {
                    // InternalElkGraphJson.g:1348:2: ( 'null' )
                    // InternalElkGraphJson.g:1349:3: 'null'
                    {
                     before(grammarAccess.getJsonValueAccess().getNullKeyword_5()); 
                    match(input,15,FOLLOW_2); 
                     after(grammarAccess.getJsonValueAccess().getNullKeyword_5()); 

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
    // $ANTLR end "rule__JsonValue__Alternatives"


    // $ANTLR start "rule__KeyChildren__Alternatives"
    // InternalElkGraphJson.g:1358:1: rule__KeyChildren__Alternatives : ( ( '\"children\"' ) | ( '\\'children\\'' ) | ( 'children' ) );
    public final void rule__KeyChildren__Alternatives() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:1362:1: ( ( '\"children\"' ) | ( '\\'children\\'' ) | ( 'children' ) )
            int alt13=3;
            switch ( input.LA(1) ) {
            case 16:
                {
                alt13=1;
                }
                break;
            case 17:
                {
                alt13=2;
                }
                break;
            case 18:
                {
                alt13=3;
                }
                break;
            default:
                NoViableAltException nvae =
                    new NoViableAltException("", 13, 0, input);

                throw nvae;
            }

            switch (alt13) {
                case 1 :
                    // InternalElkGraphJson.g:1363:2: ( '\"children\"' )
                    {
                    // InternalElkGraphJson.g:1363:2: ( '\"children\"' )
                    // InternalElkGraphJson.g:1364:3: '\"children\"'
                    {
                     before(grammarAccess.getKeyChildrenAccess().getChildrenKeyword_0()); 
                    match(input,16,FOLLOW_2); 
                     after(grammarAccess.getKeyChildrenAccess().getChildrenKeyword_0()); 

                    }


                    }
                    break;
                case 2 :
                    // InternalElkGraphJson.g:1369:2: ( '\\'children\\'' )
                    {
                    // InternalElkGraphJson.g:1369:2: ( '\\'children\\'' )
                    // InternalElkGraphJson.g:1370:3: '\\'children\\''
                    {
                     before(grammarAccess.getKeyChildrenAccess().getChildrenKeyword_1()); 
                    match(input,17,FOLLOW_2); 
                     after(grammarAccess.getKeyChildrenAccess().getChildrenKeyword_1()); 

                    }


                    }
                    break;
                case 3 :
                    // InternalElkGraphJson.g:1375:2: ( 'children' )
                    {
                    // InternalElkGraphJson.g:1375:2: ( 'children' )
                    // InternalElkGraphJson.g:1376:3: 'children'
                    {
                     before(grammarAccess.getKeyChildrenAccess().getChildrenKeyword_2()); 
                    match(input,18,FOLLOW_2); 
                     after(grammarAccess.getKeyChildrenAccess().getChildrenKeyword_2()); 

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
    // $ANTLR end "rule__KeyChildren__Alternatives"


    // $ANTLR start "rule__KeyPorts__Alternatives"
    // InternalElkGraphJson.g:1385:1: rule__KeyPorts__Alternatives : ( ( '\"ports\"' ) | ( '\\'ports\\'' ) | ( 'ports' ) );
    public final void rule__KeyPorts__Alternatives() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:1389:1: ( ( '\"ports\"' ) | ( '\\'ports\\'' ) | ( 'ports' ) )
            int alt14=3;
            switch ( input.LA(1) ) {
            case 19:
                {
                alt14=1;
                }
                break;
            case 20:
                {
                alt14=2;
                }
                break;
            case 21:
                {
                alt14=3;
                }
                break;
            default:
                NoViableAltException nvae =
                    new NoViableAltException("", 14, 0, input);

                throw nvae;
            }

            switch (alt14) {
                case 1 :
                    // InternalElkGraphJson.g:1390:2: ( '\"ports\"' )
                    {
                    // InternalElkGraphJson.g:1390:2: ( '\"ports\"' )
                    // InternalElkGraphJson.g:1391:3: '\"ports\"'
                    {
                     before(grammarAccess.getKeyPortsAccess().getPortsKeyword_0()); 
                    match(input,19,FOLLOW_2); 
                     after(grammarAccess.getKeyPortsAccess().getPortsKeyword_0()); 

                    }


                    }
                    break;
                case 2 :
                    // InternalElkGraphJson.g:1396:2: ( '\\'ports\\'' )
                    {
                    // InternalElkGraphJson.g:1396:2: ( '\\'ports\\'' )
                    // InternalElkGraphJson.g:1397:3: '\\'ports\\''
                    {
                     before(grammarAccess.getKeyPortsAccess().getPortsKeyword_1()); 
                    match(input,20,FOLLOW_2); 
                     after(grammarAccess.getKeyPortsAccess().getPortsKeyword_1()); 

                    }


                    }
                    break;
                case 3 :
                    // InternalElkGraphJson.g:1402:2: ( 'ports' )
                    {
                    // InternalElkGraphJson.g:1402:2: ( 'ports' )
                    // InternalElkGraphJson.g:1403:3: 'ports'
                    {
                     before(grammarAccess.getKeyPortsAccess().getPortsKeyword_2()); 
                    match(input,21,FOLLOW_2); 
                     after(grammarAccess.getKeyPortsAccess().getPortsKeyword_2()); 

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
    // $ANTLR end "rule__KeyPorts__Alternatives"


    // $ANTLR start "rule__KeyLabels__Alternatives"
    // InternalElkGraphJson.g:1412:1: rule__KeyLabels__Alternatives : ( ( '\"labels\"' ) | ( '\\'labels\\'' ) | ( 'labels' ) );
    public final void rule__KeyLabels__Alternatives() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:1416:1: ( ( '\"labels\"' ) | ( '\\'labels\\'' ) | ( 'labels' ) )
            int alt15=3;
            switch ( input.LA(1) ) {
            case 22:
                {
                alt15=1;
                }
                break;
            case 23:
                {
                alt15=2;
                }
                break;
            case 24:
                {
                alt15=3;
                }
                break;
            default:
                NoViableAltException nvae =
                    new NoViableAltException("", 15, 0, input);

                throw nvae;
            }

            switch (alt15) {
                case 1 :
                    // InternalElkGraphJson.g:1417:2: ( '\"labels\"' )
                    {
                    // InternalElkGraphJson.g:1417:2: ( '\"labels\"' )
                    // InternalElkGraphJson.g:1418:3: '\"labels\"'
                    {
                     before(grammarAccess.getKeyLabelsAccess().getLabelsKeyword_0()); 
                    match(input,22,FOLLOW_2); 
                     after(grammarAccess.getKeyLabelsAccess().getLabelsKeyword_0()); 

                    }


                    }
                    break;
                case 2 :
                    // InternalElkGraphJson.g:1423:2: ( '\\'labels\\'' )
                    {
                    // InternalElkGraphJson.g:1423:2: ( '\\'labels\\'' )
                    // InternalElkGraphJson.g:1424:3: '\\'labels\\''
                    {
                     before(grammarAccess.getKeyLabelsAccess().getLabelsKeyword_1()); 
                    match(input,23,FOLLOW_2); 
                     after(grammarAccess.getKeyLabelsAccess().getLabelsKeyword_1()); 

                    }


                    }
                    break;
                case 3 :
                    // InternalElkGraphJson.g:1429:2: ( 'labels' )
                    {
                    // InternalElkGraphJson.g:1429:2: ( 'labels' )
                    // InternalElkGraphJson.g:1430:3: 'labels'
                    {
                     before(grammarAccess.getKeyLabelsAccess().getLabelsKeyword_2()); 
                    match(input,24,FOLLOW_2); 
                     after(grammarAccess.getKeyLabelsAccess().getLabelsKeyword_2()); 

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
    // $ANTLR end "rule__KeyLabels__Alternatives"


    // $ANTLR start "rule__KeyEdges__Alternatives"
    // InternalElkGraphJson.g:1439:1: rule__KeyEdges__Alternatives : ( ( '\"edges\"' ) | ( '\\'edges\\'' ) | ( 'edges' ) );
    public final void rule__KeyEdges__Alternatives() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:1443:1: ( ( '\"edges\"' ) | ( '\\'edges\\'' ) | ( 'edges' ) )
            int alt16=3;
            switch ( input.LA(1) ) {
            case 25:
                {
                alt16=1;
                }
                break;
            case 26:
                {
                alt16=2;
                }
                break;
            case 27:
                {
                alt16=3;
                }
                break;
            default:
                NoViableAltException nvae =
                    new NoViableAltException("", 16, 0, input);

                throw nvae;
            }

            switch (alt16) {
                case 1 :
                    // InternalElkGraphJson.g:1444:2: ( '\"edges\"' )
                    {
                    // InternalElkGraphJson.g:1444:2: ( '\"edges\"' )
                    // InternalElkGraphJson.g:1445:3: '\"edges\"'
                    {
                     before(grammarAccess.getKeyEdgesAccess().getEdgesKeyword_0()); 
                    match(input,25,FOLLOW_2); 
                     after(grammarAccess.getKeyEdgesAccess().getEdgesKeyword_0()); 

                    }


                    }
                    break;
                case 2 :
                    // InternalElkGraphJson.g:1450:2: ( '\\'edges\\'' )
                    {
                    // InternalElkGraphJson.g:1450:2: ( '\\'edges\\'' )
                    // InternalElkGraphJson.g:1451:3: '\\'edges\\''
                    {
                     before(grammarAccess.getKeyEdgesAccess().getEdgesKeyword_1()); 
                    match(input,26,FOLLOW_2); 
                     after(grammarAccess.getKeyEdgesAccess().getEdgesKeyword_1()); 

                    }


                    }
                    break;
                case 3 :
                    // InternalElkGraphJson.g:1456:2: ( 'edges' )
                    {
                    // InternalElkGraphJson.g:1456:2: ( 'edges' )
                    // InternalElkGraphJson.g:1457:3: 'edges'
                    {
                     before(grammarAccess.getKeyEdgesAccess().getEdgesKeyword_2()); 
                    match(input,27,FOLLOW_2); 
                     after(grammarAccess.getKeyEdgesAccess().getEdgesKeyword_2()); 

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
    // $ANTLR end "rule__KeyEdges__Alternatives"


    // $ANTLR start "rule__KeyLayoutOptions__Alternatives"
    // InternalElkGraphJson.g:1466:1: rule__KeyLayoutOptions__Alternatives : ( ( '\"layoutOptions\"' ) | ( '\\'layoutOptions\\'' ) | ( 'layoutOptions' ) | ( '\"properties\"' ) | ( '\\'properties\\'' ) | ( 'properties' ) );
    public final void rule__KeyLayoutOptions__Alternatives() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:1470:1: ( ( '\"layoutOptions\"' ) | ( '\\'layoutOptions\\'' ) | ( 'layoutOptions' ) | ( '\"properties\"' ) | ( '\\'properties\\'' ) | ( 'properties' ) )
            int alt17=6;
            switch ( input.LA(1) ) {
            case 28:
                {
                alt17=1;
                }
                break;
            case 29:
                {
                alt17=2;
                }
                break;
            case 30:
                {
                alt17=3;
                }
                break;
            case 31:
                {
                alt17=4;
                }
                break;
            case 32:
                {
                alt17=5;
                }
                break;
            case 33:
                {
                alt17=6;
                }
                break;
            default:
                NoViableAltException nvae =
                    new NoViableAltException("", 17, 0, input);

                throw nvae;
            }

            switch (alt17) {
                case 1 :
                    // InternalElkGraphJson.g:1471:2: ( '\"layoutOptions\"' )
                    {
                    // InternalElkGraphJson.g:1471:2: ( '\"layoutOptions\"' )
                    // InternalElkGraphJson.g:1472:3: '\"layoutOptions\"'
                    {
                     before(grammarAccess.getKeyLayoutOptionsAccess().getLayoutOptionsKeyword_0()); 
                    match(input,28,FOLLOW_2); 
                     after(grammarAccess.getKeyLayoutOptionsAccess().getLayoutOptionsKeyword_0()); 

                    }


                    }
                    break;
                case 2 :
                    // InternalElkGraphJson.g:1477:2: ( '\\'layoutOptions\\'' )
                    {
                    // InternalElkGraphJson.g:1477:2: ( '\\'layoutOptions\\'' )
                    // InternalElkGraphJson.g:1478:3: '\\'layoutOptions\\''
                    {
                     before(grammarAccess.getKeyLayoutOptionsAccess().getLayoutOptionsKeyword_1()); 
                    match(input,29,FOLLOW_2); 
                     after(grammarAccess.getKeyLayoutOptionsAccess().getLayoutOptionsKeyword_1()); 

                    }


                    }
                    break;
                case 3 :
                    // InternalElkGraphJson.g:1483:2: ( 'layoutOptions' )
                    {
                    // InternalElkGraphJson.g:1483:2: ( 'layoutOptions' )
                    // InternalElkGraphJson.g:1484:3: 'layoutOptions'
                    {
                     before(grammarAccess.getKeyLayoutOptionsAccess().getLayoutOptionsKeyword_2()); 
                    match(input,30,FOLLOW_2); 
                     after(grammarAccess.getKeyLayoutOptionsAccess().getLayoutOptionsKeyword_2()); 

                    }


                    }
                    break;
                case 4 :
                    // InternalElkGraphJson.g:1489:2: ( '\"properties\"' )
                    {
                    // InternalElkGraphJson.g:1489:2: ( '\"properties\"' )
                    // InternalElkGraphJson.g:1490:3: '\"properties\"'
                    {
                     before(grammarAccess.getKeyLayoutOptionsAccess().getPropertiesKeyword_3()); 
                    match(input,31,FOLLOW_2); 
                     after(grammarAccess.getKeyLayoutOptionsAccess().getPropertiesKeyword_3()); 

                    }


                    }
                    break;
                case 5 :
                    // InternalElkGraphJson.g:1495:2: ( '\\'properties\\'' )
                    {
                    // InternalElkGraphJson.g:1495:2: ( '\\'properties\\'' )
                    // InternalElkGraphJson.g:1496:3: '\\'properties\\''
                    {
                     before(grammarAccess.getKeyLayoutOptionsAccess().getPropertiesKeyword_4()); 
                    match(input,32,FOLLOW_2); 
                     after(grammarAccess.getKeyLayoutOptionsAccess().getPropertiesKeyword_4()); 

                    }


                    }
                    break;
                case 6 :
                    // InternalElkGraphJson.g:1501:2: ( 'properties' )
                    {
                    // InternalElkGraphJson.g:1501:2: ( 'properties' )
                    // InternalElkGraphJson.g:1502:3: 'properties'
                    {
                     before(grammarAccess.getKeyLayoutOptionsAccess().getPropertiesKeyword_5()); 
                    match(input,33,FOLLOW_2); 
                     after(grammarAccess.getKeyLayoutOptionsAccess().getPropertiesKeyword_5()); 

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
    // $ANTLR end "rule__KeyLayoutOptions__Alternatives"


    // $ANTLR start "rule__KeyId__Alternatives"
    // InternalElkGraphJson.g:1511:1: rule__KeyId__Alternatives : ( ( '\"id\"' ) | ( '\\'id\\'' ) | ( 'id' ) );
    public final void rule__KeyId__Alternatives() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:1515:1: ( ( '\"id\"' ) | ( '\\'id\\'' ) | ( 'id' ) )
            int alt18=3;
            switch ( input.LA(1) ) {
            case 34:
                {
                alt18=1;
                }
                break;
            case 35:
                {
                alt18=2;
                }
                break;
            case 36:
                {
                alt18=3;
                }
                break;
            default:
                NoViableAltException nvae =
                    new NoViableAltException("", 18, 0, input);

                throw nvae;
            }

            switch (alt18) {
                case 1 :
                    // InternalElkGraphJson.g:1516:2: ( '\"id\"' )
                    {
                    // InternalElkGraphJson.g:1516:2: ( '\"id\"' )
                    // InternalElkGraphJson.g:1517:3: '\"id\"'
                    {
                     before(grammarAccess.getKeyIdAccess().getIdKeyword_0()); 
                    match(input,34,FOLLOW_2); 
                     after(grammarAccess.getKeyIdAccess().getIdKeyword_0()); 

                    }


                    }
                    break;
                case 2 :
                    // InternalElkGraphJson.g:1522:2: ( '\\'id\\'' )
                    {
                    // InternalElkGraphJson.g:1522:2: ( '\\'id\\'' )
                    // InternalElkGraphJson.g:1523:3: '\\'id\\''
                    {
                     before(grammarAccess.getKeyIdAccess().getIdKeyword_1()); 
                    match(input,35,FOLLOW_2); 
                     after(grammarAccess.getKeyIdAccess().getIdKeyword_1()); 

                    }


                    }
                    break;
                case 3 :
                    // InternalElkGraphJson.g:1528:2: ( 'id' )
                    {
                    // InternalElkGraphJson.g:1528:2: ( 'id' )
                    // InternalElkGraphJson.g:1529:3: 'id'
                    {
                     before(grammarAccess.getKeyIdAccess().getIdKeyword_2()); 
                    match(input,36,FOLLOW_2); 
                     after(grammarAccess.getKeyIdAccess().getIdKeyword_2()); 

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
    // $ANTLR end "rule__KeyId__Alternatives"


    // $ANTLR start "rule__KeyX__Alternatives"
    // InternalElkGraphJson.g:1538:1: rule__KeyX__Alternatives : ( ( '\"x\"' ) | ( '\\'x\\'' ) | ( 'x' ) );
    public final void rule__KeyX__Alternatives() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:1542:1: ( ( '\"x\"' ) | ( '\\'x\\'' ) | ( 'x' ) )
            int alt19=3;
            switch ( input.LA(1) ) {
            case 37:
                {
                alt19=1;
                }
                break;
            case 38:
                {
                alt19=2;
                }
                break;
            case 39:
                {
                alt19=3;
                }
                break;
            default:
                NoViableAltException nvae =
                    new NoViableAltException("", 19, 0, input);

                throw nvae;
            }

            switch (alt19) {
                case 1 :
                    // InternalElkGraphJson.g:1543:2: ( '\"x\"' )
                    {
                    // InternalElkGraphJson.g:1543:2: ( '\"x\"' )
                    // InternalElkGraphJson.g:1544:3: '\"x\"'
                    {
                     before(grammarAccess.getKeyXAccess().getXKeyword_0()); 
                    match(input,37,FOLLOW_2); 
                     after(grammarAccess.getKeyXAccess().getXKeyword_0()); 

                    }


                    }
                    break;
                case 2 :
                    // InternalElkGraphJson.g:1549:2: ( '\\'x\\'' )
                    {
                    // InternalElkGraphJson.g:1549:2: ( '\\'x\\'' )
                    // InternalElkGraphJson.g:1550:3: '\\'x\\''
                    {
                     before(grammarAccess.getKeyXAccess().getXKeyword_1()); 
                    match(input,38,FOLLOW_2); 
                     after(grammarAccess.getKeyXAccess().getXKeyword_1()); 

                    }


                    }
                    break;
                case 3 :
                    // InternalElkGraphJson.g:1555:2: ( 'x' )
                    {
                    // InternalElkGraphJson.g:1555:2: ( 'x' )
                    // InternalElkGraphJson.g:1556:3: 'x'
                    {
                     before(grammarAccess.getKeyXAccess().getXKeyword_2()); 
                    match(input,39,FOLLOW_2); 
                     after(grammarAccess.getKeyXAccess().getXKeyword_2()); 

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
    // $ANTLR end "rule__KeyX__Alternatives"


    // $ANTLR start "rule__KeyY__Alternatives"
    // InternalElkGraphJson.g:1565:1: rule__KeyY__Alternatives : ( ( '\"y\"' ) | ( '\\'y\\'' ) | ( 'y' ) );
    public final void rule__KeyY__Alternatives() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:1569:1: ( ( '\"y\"' ) | ( '\\'y\\'' ) | ( 'y' ) )
            int alt20=3;
            switch ( input.LA(1) ) {
            case 40:
                {
                alt20=1;
                }
                break;
            case 41:
                {
                alt20=2;
                }
                break;
            case 42:
                {
                alt20=3;
                }
                break;
            default:
                NoViableAltException nvae =
                    new NoViableAltException("", 20, 0, input);

                throw nvae;
            }

            switch (alt20) {
                case 1 :
                    // InternalElkGraphJson.g:1570:2: ( '\"y\"' )
                    {
                    // InternalElkGraphJson.g:1570:2: ( '\"y\"' )
                    // InternalElkGraphJson.g:1571:3: '\"y\"'
                    {
                     before(grammarAccess.getKeyYAccess().getYKeyword_0()); 
                    match(input,40,FOLLOW_2); 
                     after(grammarAccess.getKeyYAccess().getYKeyword_0()); 

                    }


                    }
                    break;
                case 2 :
                    // InternalElkGraphJson.g:1576:2: ( '\\'y\\'' )
                    {
                    // InternalElkGraphJson.g:1576:2: ( '\\'y\\'' )
                    // InternalElkGraphJson.g:1577:3: '\\'y\\''
                    {
                     before(grammarAccess.getKeyYAccess().getYKeyword_1()); 
                    match(input,41,FOLLOW_2); 
                     after(grammarAccess.getKeyYAccess().getYKeyword_1()); 

                    }


                    }
                    break;
                case 3 :
                    // InternalElkGraphJson.g:1582:2: ( 'y' )
                    {
                    // InternalElkGraphJson.g:1582:2: ( 'y' )
                    // InternalElkGraphJson.g:1583:3: 'y'
                    {
                     before(grammarAccess.getKeyYAccess().getYKeyword_2()); 
                    match(input,42,FOLLOW_2); 
                     after(grammarAccess.getKeyYAccess().getYKeyword_2()); 

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
    // $ANTLR end "rule__KeyY__Alternatives"


    // $ANTLR start "rule__KeyWidth__Alternatives"
    // InternalElkGraphJson.g:1592:1: rule__KeyWidth__Alternatives : ( ( '\"width\"' ) | ( '\\'width\\'' ) | ( 'width' ) );
    public final void rule__KeyWidth__Alternatives() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:1596:1: ( ( '\"width\"' ) | ( '\\'width\\'' ) | ( 'width' ) )
            int alt21=3;
            switch ( input.LA(1) ) {
            case 43:
                {
                alt21=1;
                }
                break;
            case 44:
                {
                alt21=2;
                }
                break;
            case 45:
                {
                alt21=3;
                }
                break;
            default:
                NoViableAltException nvae =
                    new NoViableAltException("", 21, 0, input);

                throw nvae;
            }

            switch (alt21) {
                case 1 :
                    // InternalElkGraphJson.g:1597:2: ( '\"width\"' )
                    {
                    // InternalElkGraphJson.g:1597:2: ( '\"width\"' )
                    // InternalElkGraphJson.g:1598:3: '\"width\"'
                    {
                     before(grammarAccess.getKeyWidthAccess().getWidthKeyword_0()); 
                    match(input,43,FOLLOW_2); 
                     after(grammarAccess.getKeyWidthAccess().getWidthKeyword_0()); 

                    }


                    }
                    break;
                case 2 :
                    // InternalElkGraphJson.g:1603:2: ( '\\'width\\'' )
                    {
                    // InternalElkGraphJson.g:1603:2: ( '\\'width\\'' )
                    // InternalElkGraphJson.g:1604:3: '\\'width\\''
                    {
                     before(grammarAccess.getKeyWidthAccess().getWidthKeyword_1()); 
                    match(input,44,FOLLOW_2); 
                     after(grammarAccess.getKeyWidthAccess().getWidthKeyword_1()); 

                    }


                    }
                    break;
                case 3 :
                    // InternalElkGraphJson.g:1609:2: ( 'width' )
                    {
                    // InternalElkGraphJson.g:1609:2: ( 'width' )
                    // InternalElkGraphJson.g:1610:3: 'width'
                    {
                     before(grammarAccess.getKeyWidthAccess().getWidthKeyword_2()); 
                    match(input,45,FOLLOW_2); 
                     after(grammarAccess.getKeyWidthAccess().getWidthKeyword_2()); 

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
    // $ANTLR end "rule__KeyWidth__Alternatives"


    // $ANTLR start "rule__KeyHeight__Alternatives"
    // InternalElkGraphJson.g:1619:1: rule__KeyHeight__Alternatives : ( ( '\"height\"' ) | ( '\\'height\\'' ) | ( 'height' ) );
    public final void rule__KeyHeight__Alternatives() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:1623:1: ( ( '\"height\"' ) | ( '\\'height\\'' ) | ( 'height' ) )
            int alt22=3;
            switch ( input.LA(1) ) {
            case 46:
                {
                alt22=1;
                }
                break;
            case 47:
                {
                alt22=2;
                }
                break;
            case 48:
                {
                alt22=3;
                }
                break;
            default:
                NoViableAltException nvae =
                    new NoViableAltException("", 22, 0, input);

                throw nvae;
            }

            switch (alt22) {
                case 1 :
                    // InternalElkGraphJson.g:1624:2: ( '\"height\"' )
                    {
                    // InternalElkGraphJson.g:1624:2: ( '\"height\"' )
                    // InternalElkGraphJson.g:1625:3: '\"height\"'
                    {
                     before(grammarAccess.getKeyHeightAccess().getHeightKeyword_0()); 
                    match(input,46,FOLLOW_2); 
                     after(grammarAccess.getKeyHeightAccess().getHeightKeyword_0()); 

                    }


                    }
                    break;
                case 2 :
                    // InternalElkGraphJson.g:1630:2: ( '\\'height\\'' )
                    {
                    // InternalElkGraphJson.g:1630:2: ( '\\'height\\'' )
                    // InternalElkGraphJson.g:1631:3: '\\'height\\''
                    {
                     before(grammarAccess.getKeyHeightAccess().getHeightKeyword_1()); 
                    match(input,47,FOLLOW_2); 
                     after(grammarAccess.getKeyHeightAccess().getHeightKeyword_1()); 

                    }


                    }
                    break;
                case 3 :
                    // InternalElkGraphJson.g:1636:2: ( 'height' )
                    {
                    // InternalElkGraphJson.g:1636:2: ( 'height' )
                    // InternalElkGraphJson.g:1637:3: 'height'
                    {
                     before(grammarAccess.getKeyHeightAccess().getHeightKeyword_2()); 
                    match(input,48,FOLLOW_2); 
                     after(grammarAccess.getKeyHeightAccess().getHeightKeyword_2()); 

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
    // $ANTLR end "rule__KeyHeight__Alternatives"


    // $ANTLR start "rule__KeySources__Alternatives"
    // InternalElkGraphJson.g:1646:1: rule__KeySources__Alternatives : ( ( '\"sources\"' ) | ( '\\'sources\\'' ) | ( 'sources' ) );
    public final void rule__KeySources__Alternatives() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:1650:1: ( ( '\"sources\"' ) | ( '\\'sources\\'' ) | ( 'sources' ) )
            int alt23=3;
            switch ( input.LA(1) ) {
            case 49:
                {
                alt23=1;
                }
                break;
            case 50:
                {
                alt23=2;
                }
                break;
            case 51:
                {
                alt23=3;
                }
                break;
            default:
                NoViableAltException nvae =
                    new NoViableAltException("", 23, 0, input);

                throw nvae;
            }

            switch (alt23) {
                case 1 :
                    // InternalElkGraphJson.g:1651:2: ( '\"sources\"' )
                    {
                    // InternalElkGraphJson.g:1651:2: ( '\"sources\"' )
                    // InternalElkGraphJson.g:1652:3: '\"sources\"'
                    {
                     before(grammarAccess.getKeySourcesAccess().getSourcesKeyword_0()); 
                    match(input,49,FOLLOW_2); 
                     after(grammarAccess.getKeySourcesAccess().getSourcesKeyword_0()); 

                    }


                    }
                    break;
                case 2 :
                    // InternalElkGraphJson.g:1657:2: ( '\\'sources\\'' )
                    {
                    // InternalElkGraphJson.g:1657:2: ( '\\'sources\\'' )
                    // InternalElkGraphJson.g:1658:3: '\\'sources\\''
                    {
                     before(grammarAccess.getKeySourcesAccess().getSourcesKeyword_1()); 
                    match(input,50,FOLLOW_2); 
                     after(grammarAccess.getKeySourcesAccess().getSourcesKeyword_1()); 

                    }


                    }
                    break;
                case 3 :
                    // InternalElkGraphJson.g:1663:2: ( 'sources' )
                    {
                    // InternalElkGraphJson.g:1663:2: ( 'sources' )
                    // InternalElkGraphJson.g:1664:3: 'sources'
                    {
                     before(grammarAccess.getKeySourcesAccess().getSourcesKeyword_2()); 
                    match(input,51,FOLLOW_2); 
                     after(grammarAccess.getKeySourcesAccess().getSourcesKeyword_2()); 

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
    // $ANTLR end "rule__KeySources__Alternatives"


    // $ANTLR start "rule__KeyTargets__Alternatives"
    // InternalElkGraphJson.g:1673:1: rule__KeyTargets__Alternatives : ( ( '\"targets\"' ) | ( '\\'targets\\'' ) | ( 'targets' ) );
    public final void rule__KeyTargets__Alternatives() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:1677:1: ( ( '\"targets\"' ) | ( '\\'targets\\'' ) | ( 'targets' ) )
            int alt24=3;
            switch ( input.LA(1) ) {
            case 52:
                {
                alt24=1;
                }
                break;
            case 53:
                {
                alt24=2;
                }
                break;
            case 54:
                {
                alt24=3;
                }
                break;
            default:
                NoViableAltException nvae =
                    new NoViableAltException("", 24, 0, input);

                throw nvae;
            }

            switch (alt24) {
                case 1 :
                    // InternalElkGraphJson.g:1678:2: ( '\"targets\"' )
                    {
                    // InternalElkGraphJson.g:1678:2: ( '\"targets\"' )
                    // InternalElkGraphJson.g:1679:3: '\"targets\"'
                    {
                     before(grammarAccess.getKeyTargetsAccess().getTargetsKeyword_0()); 
                    match(input,52,FOLLOW_2); 
                     after(grammarAccess.getKeyTargetsAccess().getTargetsKeyword_0()); 

                    }


                    }
                    break;
                case 2 :
                    // InternalElkGraphJson.g:1684:2: ( '\\'targets\\'' )
                    {
                    // InternalElkGraphJson.g:1684:2: ( '\\'targets\\'' )
                    // InternalElkGraphJson.g:1685:3: '\\'targets\\''
                    {
                     before(grammarAccess.getKeyTargetsAccess().getTargetsKeyword_1()); 
                    match(input,53,FOLLOW_2); 
                     after(grammarAccess.getKeyTargetsAccess().getTargetsKeyword_1()); 

                    }


                    }
                    break;
                case 3 :
                    // InternalElkGraphJson.g:1690:2: ( 'targets' )
                    {
                    // InternalElkGraphJson.g:1690:2: ( 'targets' )
                    // InternalElkGraphJson.g:1691:3: 'targets'
                    {
                     before(grammarAccess.getKeyTargetsAccess().getTargetsKeyword_2()); 
                    match(input,54,FOLLOW_2); 
                     after(grammarAccess.getKeyTargetsAccess().getTargetsKeyword_2()); 

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
    // $ANTLR end "rule__KeyTargets__Alternatives"


    // $ANTLR start "rule__KeyText__Alternatives"
    // InternalElkGraphJson.g:1700:1: rule__KeyText__Alternatives : ( ( '\"text\"' ) | ( '\\'text\\'' ) | ( 'text' ) );
    public final void rule__KeyText__Alternatives() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:1704:1: ( ( '\"text\"' ) | ( '\\'text\\'' ) | ( 'text' ) )
            int alt25=3;
            switch ( input.LA(1) ) {
            case 55:
                {
                alt25=1;
                }
                break;
            case 56:
                {
                alt25=2;
                }
                break;
            case 57:
                {
                alt25=3;
                }
                break;
            default:
                NoViableAltException nvae =
                    new NoViableAltException("", 25, 0, input);

                throw nvae;
            }

            switch (alt25) {
                case 1 :
                    // InternalElkGraphJson.g:1705:2: ( '\"text\"' )
                    {
                    // InternalElkGraphJson.g:1705:2: ( '\"text\"' )
                    // InternalElkGraphJson.g:1706:3: '\"text\"'
                    {
                     before(grammarAccess.getKeyTextAccess().getTextKeyword_0()); 
                    match(input,55,FOLLOW_2); 
                     after(grammarAccess.getKeyTextAccess().getTextKeyword_0()); 

                    }


                    }
                    break;
                case 2 :
                    // InternalElkGraphJson.g:1711:2: ( '\\'text\\'' )
                    {
                    // InternalElkGraphJson.g:1711:2: ( '\\'text\\'' )
                    // InternalElkGraphJson.g:1712:3: '\\'text\\''
                    {
                     before(grammarAccess.getKeyTextAccess().getTextKeyword_1()); 
                    match(input,56,FOLLOW_2); 
                     after(grammarAccess.getKeyTextAccess().getTextKeyword_1()); 

                    }


                    }
                    break;
                case 3 :
                    // InternalElkGraphJson.g:1717:2: ( 'text' )
                    {
                    // InternalElkGraphJson.g:1717:2: ( 'text' )
                    // InternalElkGraphJson.g:1718:3: 'text'
                    {
                     before(grammarAccess.getKeyTextAccess().getTextKeyword_2()); 
                    match(input,57,FOLLOW_2); 
                     after(grammarAccess.getKeyTextAccess().getTextKeyword_2()); 

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
    // $ANTLR end "rule__KeyText__Alternatives"


    // $ANTLR start "rule__ElkNode__Group__0"
    // InternalElkGraphJson.g:1727:1: rule__ElkNode__Group__0 : rule__ElkNode__Group__0__Impl rule__ElkNode__Group__1 ;
    public final void rule__ElkNode__Group__0() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:1731:1: ( rule__ElkNode__Group__0__Impl rule__ElkNode__Group__1 )
            // InternalElkGraphJson.g:1732:2: rule__ElkNode__Group__0__Impl rule__ElkNode__Group__1
            {
            pushFollow(FOLLOW_3);
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
    // InternalElkGraphJson.g:1739:1: rule__ElkNode__Group__0__Impl : ( () ) ;
    public final void rule__ElkNode__Group__0__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:1743:1: ( ( () ) )
            // InternalElkGraphJson.g:1744:1: ( () )
            {
            // InternalElkGraphJson.g:1744:1: ( () )
            // InternalElkGraphJson.g:1745:2: ()
            {
             before(grammarAccess.getElkNodeAccess().getElkNodeAction_0()); 
            // InternalElkGraphJson.g:1746:2: ()
            // InternalElkGraphJson.g:1746:3: 
            {
            }

             after(grammarAccess.getElkNodeAccess().getElkNodeAction_0()); 

            }


            }

        }
        finally {

            	restoreStackSize(stackSize);

        }
        return ;
    }
    // $ANTLR end "rule__ElkNode__Group__0__Impl"


    // $ANTLR start "rule__ElkNode__Group__1"
    // InternalElkGraphJson.g:1754:1: rule__ElkNode__Group__1 : rule__ElkNode__Group__1__Impl rule__ElkNode__Group__2 ;
    public final void rule__ElkNode__Group__1() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:1758:1: ( rule__ElkNode__Group__1__Impl rule__ElkNode__Group__2 )
            // InternalElkGraphJson.g:1759:2: rule__ElkNode__Group__1__Impl rule__ElkNode__Group__2
            {
            pushFollow(FOLLOW_4);
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
    // InternalElkGraphJson.g:1766:1: rule__ElkNode__Group__1__Impl : ( '{' ) ;
    public final void rule__ElkNode__Group__1__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:1770:1: ( ( '{' ) )
            // InternalElkGraphJson.g:1771:1: ( '{' )
            {
            // InternalElkGraphJson.g:1771:1: ( '{' )
            // InternalElkGraphJson.g:1772:2: '{'
            {
             before(grammarAccess.getElkNodeAccess().getLeftCurlyBracketKeyword_1()); 
            match(input,58,FOLLOW_2); 
             after(grammarAccess.getElkNodeAccess().getLeftCurlyBracketKeyword_1()); 

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
    // InternalElkGraphJson.g:1781:1: rule__ElkNode__Group__2 : rule__ElkNode__Group__2__Impl rule__ElkNode__Group__3 ;
    public final void rule__ElkNode__Group__2() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:1785:1: ( rule__ElkNode__Group__2__Impl rule__ElkNode__Group__3 )
            // InternalElkGraphJson.g:1786:2: rule__ElkNode__Group__2__Impl rule__ElkNode__Group__3
            {
            pushFollow(FOLLOW_4);
            rule__ElkNode__Group__2__Impl();

            state._fsp--;

            pushFollow(FOLLOW_2);
            rule__ElkNode__Group__3();

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
    // InternalElkGraphJson.g:1793:1: rule__ElkNode__Group__2__Impl : ( ( rule__ElkNode__Group_2__0 )? ) ;
    public final void rule__ElkNode__Group__2__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:1797:1: ( ( ( rule__ElkNode__Group_2__0 )? ) )
            // InternalElkGraphJson.g:1798:1: ( ( rule__ElkNode__Group_2__0 )? )
            {
            // InternalElkGraphJson.g:1798:1: ( ( rule__ElkNode__Group_2__0 )? )
            // InternalElkGraphJson.g:1799:2: ( rule__ElkNode__Group_2__0 )?
            {
             before(grammarAccess.getElkNodeAccess().getGroup_2()); 
            // InternalElkGraphJson.g:1800:2: ( rule__ElkNode__Group_2__0 )?
            int alt26=2;
            int LA26_0 = input.LA(1);

            if ( ((LA26_0>=RULE_STRING && LA26_0<=RULE_ID)||(LA26_0>=16 && LA26_0<=48)) ) {
                alt26=1;
            }
            switch (alt26) {
                case 1 :
                    // InternalElkGraphJson.g:1800:3: rule__ElkNode__Group_2__0
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


    // $ANTLR start "rule__ElkNode__Group__3"
    // InternalElkGraphJson.g:1808:1: rule__ElkNode__Group__3 : rule__ElkNode__Group__3__Impl rule__ElkNode__Group__4 ;
    public final void rule__ElkNode__Group__3() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:1812:1: ( rule__ElkNode__Group__3__Impl rule__ElkNode__Group__4 )
            // InternalElkGraphJson.g:1813:2: rule__ElkNode__Group__3__Impl rule__ElkNode__Group__4
            {
            pushFollow(FOLLOW_4);
            rule__ElkNode__Group__3__Impl();

            state._fsp--;

            pushFollow(FOLLOW_2);
            rule__ElkNode__Group__4();

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
    // $ANTLR end "rule__ElkNode__Group__3"


    // $ANTLR start "rule__ElkNode__Group__3__Impl"
    // InternalElkGraphJson.g:1820:1: rule__ElkNode__Group__3__Impl : ( ( ',' )? ) ;
    public final void rule__ElkNode__Group__3__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:1824:1: ( ( ( ',' )? ) )
            // InternalElkGraphJson.g:1825:1: ( ( ',' )? )
            {
            // InternalElkGraphJson.g:1825:1: ( ( ',' )? )
            // InternalElkGraphJson.g:1826:2: ( ',' )?
            {
             before(grammarAccess.getElkNodeAccess().getCommaKeyword_3()); 
            // InternalElkGraphJson.g:1827:2: ( ',' )?
            int alt27=2;
            int LA27_0 = input.LA(1);

            if ( (LA27_0==59) ) {
                alt27=1;
            }
            switch (alt27) {
                case 1 :
                    // InternalElkGraphJson.g:1827:3: ','
                    {
                    match(input,59,FOLLOW_2); 

                    }
                    break;

            }

             after(grammarAccess.getElkNodeAccess().getCommaKeyword_3()); 

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
    // $ANTLR end "rule__ElkNode__Group__3__Impl"


    // $ANTLR start "rule__ElkNode__Group__4"
    // InternalElkGraphJson.g:1835:1: rule__ElkNode__Group__4 : rule__ElkNode__Group__4__Impl ;
    public final void rule__ElkNode__Group__4() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:1839:1: ( rule__ElkNode__Group__4__Impl )
            // InternalElkGraphJson.g:1840:2: rule__ElkNode__Group__4__Impl
            {
            pushFollow(FOLLOW_2);
            rule__ElkNode__Group__4__Impl();

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
    // $ANTLR end "rule__ElkNode__Group__4"


    // $ANTLR start "rule__ElkNode__Group__4__Impl"
    // InternalElkGraphJson.g:1846:1: rule__ElkNode__Group__4__Impl : ( '}' ) ;
    public final void rule__ElkNode__Group__4__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:1850:1: ( ( '}' ) )
            // InternalElkGraphJson.g:1851:1: ( '}' )
            {
            // InternalElkGraphJson.g:1851:1: ( '}' )
            // InternalElkGraphJson.g:1852:2: '}'
            {
             before(grammarAccess.getElkNodeAccess().getRightCurlyBracketKeyword_4()); 
            match(input,60,FOLLOW_2); 
             after(grammarAccess.getElkNodeAccess().getRightCurlyBracketKeyword_4()); 

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
    // $ANTLR end "rule__ElkNode__Group__4__Impl"


    // $ANTLR start "rule__ElkNode__Group_2__0"
    // InternalElkGraphJson.g:1862:1: rule__ElkNode__Group_2__0 : rule__ElkNode__Group_2__0__Impl rule__ElkNode__Group_2__1 ;
    public final void rule__ElkNode__Group_2__0() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:1866:1: ( rule__ElkNode__Group_2__0__Impl rule__ElkNode__Group_2__1 )
            // InternalElkGraphJson.g:1867:2: rule__ElkNode__Group_2__0__Impl rule__ElkNode__Group_2__1
            {
            pushFollow(FOLLOW_5);
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
    // InternalElkGraphJson.g:1874:1: rule__ElkNode__Group_2__0__Impl : ( ruleNodeElement ) ;
    public final void rule__ElkNode__Group_2__0__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:1878:1: ( ( ruleNodeElement ) )
            // InternalElkGraphJson.g:1879:1: ( ruleNodeElement )
            {
            // InternalElkGraphJson.g:1879:1: ( ruleNodeElement )
            // InternalElkGraphJson.g:1880:2: ruleNodeElement
            {
             before(grammarAccess.getElkNodeAccess().getNodeElementParserRuleCall_2_0()); 
            pushFollow(FOLLOW_2);
            ruleNodeElement();

            state._fsp--;

             after(grammarAccess.getElkNodeAccess().getNodeElementParserRuleCall_2_0()); 

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
    // InternalElkGraphJson.g:1889:1: rule__ElkNode__Group_2__1 : rule__ElkNode__Group_2__1__Impl ;
    public final void rule__ElkNode__Group_2__1() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:1893:1: ( rule__ElkNode__Group_2__1__Impl )
            // InternalElkGraphJson.g:1894:2: rule__ElkNode__Group_2__1__Impl
            {
            pushFollow(FOLLOW_2);
            rule__ElkNode__Group_2__1__Impl();

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
    // InternalElkGraphJson.g:1900:1: rule__ElkNode__Group_2__1__Impl : ( ( rule__ElkNode__Group_2_1__0 )* ) ;
    public final void rule__ElkNode__Group_2__1__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:1904:1: ( ( ( rule__ElkNode__Group_2_1__0 )* ) )
            // InternalElkGraphJson.g:1905:1: ( ( rule__ElkNode__Group_2_1__0 )* )
            {
            // InternalElkGraphJson.g:1905:1: ( ( rule__ElkNode__Group_2_1__0 )* )
            // InternalElkGraphJson.g:1906:2: ( rule__ElkNode__Group_2_1__0 )*
            {
             before(grammarAccess.getElkNodeAccess().getGroup_2_1()); 
            // InternalElkGraphJson.g:1907:2: ( rule__ElkNode__Group_2_1__0 )*
            loop28:
            do {
                int alt28=2;
                int LA28_0 = input.LA(1);

                if ( (LA28_0==59) ) {
                    int LA28_1 = input.LA(2);

                    if ( ((LA28_1>=RULE_STRING && LA28_1<=RULE_ID)||(LA28_1>=16 && LA28_1<=48)) ) {
                        alt28=1;
                    }


                }


                switch (alt28) {
            	case 1 :
            	    // InternalElkGraphJson.g:1907:3: rule__ElkNode__Group_2_1__0
            	    {
            	    pushFollow(FOLLOW_6);
            	    rule__ElkNode__Group_2_1__0();

            	    state._fsp--;


            	    }
            	    break;

            	default :
            	    break loop28;
                }
            } while (true);

             after(grammarAccess.getElkNodeAccess().getGroup_2_1()); 

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


    // $ANTLR start "rule__ElkNode__Group_2_1__0"
    // InternalElkGraphJson.g:1916:1: rule__ElkNode__Group_2_1__0 : rule__ElkNode__Group_2_1__0__Impl rule__ElkNode__Group_2_1__1 ;
    public final void rule__ElkNode__Group_2_1__0() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:1920:1: ( rule__ElkNode__Group_2_1__0__Impl rule__ElkNode__Group_2_1__1 )
            // InternalElkGraphJson.g:1921:2: rule__ElkNode__Group_2_1__0__Impl rule__ElkNode__Group_2_1__1
            {
            pushFollow(FOLLOW_7);
            rule__ElkNode__Group_2_1__0__Impl();

            state._fsp--;

            pushFollow(FOLLOW_2);
            rule__ElkNode__Group_2_1__1();

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
    // $ANTLR end "rule__ElkNode__Group_2_1__0"


    // $ANTLR start "rule__ElkNode__Group_2_1__0__Impl"
    // InternalElkGraphJson.g:1928:1: rule__ElkNode__Group_2_1__0__Impl : ( ',' ) ;
    public final void rule__ElkNode__Group_2_1__0__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:1932:1: ( ( ',' ) )
            // InternalElkGraphJson.g:1933:1: ( ',' )
            {
            // InternalElkGraphJson.g:1933:1: ( ',' )
            // InternalElkGraphJson.g:1934:2: ','
            {
             before(grammarAccess.getElkNodeAccess().getCommaKeyword_2_1_0()); 
            match(input,59,FOLLOW_2); 
             after(grammarAccess.getElkNodeAccess().getCommaKeyword_2_1_0()); 

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
    // $ANTLR end "rule__ElkNode__Group_2_1__0__Impl"


    // $ANTLR start "rule__ElkNode__Group_2_1__1"
    // InternalElkGraphJson.g:1943:1: rule__ElkNode__Group_2_1__1 : rule__ElkNode__Group_2_1__1__Impl ;
    public final void rule__ElkNode__Group_2_1__1() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:1947:1: ( rule__ElkNode__Group_2_1__1__Impl )
            // InternalElkGraphJson.g:1948:2: rule__ElkNode__Group_2_1__1__Impl
            {
            pushFollow(FOLLOW_2);
            rule__ElkNode__Group_2_1__1__Impl();

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
    // $ANTLR end "rule__ElkNode__Group_2_1__1"


    // $ANTLR start "rule__ElkNode__Group_2_1__1__Impl"
    // InternalElkGraphJson.g:1954:1: rule__ElkNode__Group_2_1__1__Impl : ( ruleNodeElement ) ;
    public final void rule__ElkNode__Group_2_1__1__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:1958:1: ( ( ruleNodeElement ) )
            // InternalElkGraphJson.g:1959:1: ( ruleNodeElement )
            {
            // InternalElkGraphJson.g:1959:1: ( ruleNodeElement )
            // InternalElkGraphJson.g:1960:2: ruleNodeElement
            {
             before(grammarAccess.getElkNodeAccess().getNodeElementParserRuleCall_2_1_1()); 
            pushFollow(FOLLOW_2);
            ruleNodeElement();

            state._fsp--;

             after(grammarAccess.getElkNodeAccess().getNodeElementParserRuleCall_2_1_1()); 

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
    // $ANTLR end "rule__ElkNode__Group_2_1__1__Impl"


    // $ANTLR start "rule__NodeElement__Group_2__0"
    // InternalElkGraphJson.g:1970:1: rule__NodeElement__Group_2__0 : rule__NodeElement__Group_2__0__Impl rule__NodeElement__Group_2__1 ;
    public final void rule__NodeElement__Group_2__0() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:1974:1: ( rule__NodeElement__Group_2__0__Impl rule__NodeElement__Group_2__1 )
            // InternalElkGraphJson.g:1975:2: rule__NodeElement__Group_2__0__Impl rule__NodeElement__Group_2__1
            {
            pushFollow(FOLLOW_8);
            rule__NodeElement__Group_2__0__Impl();

            state._fsp--;

            pushFollow(FOLLOW_2);
            rule__NodeElement__Group_2__1();

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
    // $ANTLR end "rule__NodeElement__Group_2__0"


    // $ANTLR start "rule__NodeElement__Group_2__0__Impl"
    // InternalElkGraphJson.g:1982:1: rule__NodeElement__Group_2__0__Impl : ( ruleKeyChildren ) ;
    public final void rule__NodeElement__Group_2__0__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:1986:1: ( ( ruleKeyChildren ) )
            // InternalElkGraphJson.g:1987:1: ( ruleKeyChildren )
            {
            // InternalElkGraphJson.g:1987:1: ( ruleKeyChildren )
            // InternalElkGraphJson.g:1988:2: ruleKeyChildren
            {
             before(grammarAccess.getNodeElementAccess().getKeyChildrenParserRuleCall_2_0()); 
            pushFollow(FOLLOW_2);
            ruleKeyChildren();

            state._fsp--;

             after(grammarAccess.getNodeElementAccess().getKeyChildrenParserRuleCall_2_0()); 

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
    // $ANTLR end "rule__NodeElement__Group_2__0__Impl"


    // $ANTLR start "rule__NodeElement__Group_2__1"
    // InternalElkGraphJson.g:1997:1: rule__NodeElement__Group_2__1 : rule__NodeElement__Group_2__1__Impl rule__NodeElement__Group_2__2 ;
    public final void rule__NodeElement__Group_2__1() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:2001:1: ( rule__NodeElement__Group_2__1__Impl rule__NodeElement__Group_2__2 )
            // InternalElkGraphJson.g:2002:2: rule__NodeElement__Group_2__1__Impl rule__NodeElement__Group_2__2
            {
            pushFollow(FOLLOW_9);
            rule__NodeElement__Group_2__1__Impl();

            state._fsp--;

            pushFollow(FOLLOW_2);
            rule__NodeElement__Group_2__2();

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
    // $ANTLR end "rule__NodeElement__Group_2__1"


    // $ANTLR start "rule__NodeElement__Group_2__1__Impl"
    // InternalElkGraphJson.g:2009:1: rule__NodeElement__Group_2__1__Impl : ( ':' ) ;
    public final void rule__NodeElement__Group_2__1__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:2013:1: ( ( ':' ) )
            // InternalElkGraphJson.g:2014:1: ( ':' )
            {
            // InternalElkGraphJson.g:2014:1: ( ':' )
            // InternalElkGraphJson.g:2015:2: ':'
            {
             before(grammarAccess.getNodeElementAccess().getColonKeyword_2_1()); 
            match(input,61,FOLLOW_2); 
             after(grammarAccess.getNodeElementAccess().getColonKeyword_2_1()); 

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
    // $ANTLR end "rule__NodeElement__Group_2__1__Impl"


    // $ANTLR start "rule__NodeElement__Group_2__2"
    // InternalElkGraphJson.g:2024:1: rule__NodeElement__Group_2__2 : rule__NodeElement__Group_2__2__Impl ;
    public final void rule__NodeElement__Group_2__2() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:2028:1: ( rule__NodeElement__Group_2__2__Impl )
            // InternalElkGraphJson.g:2029:2: rule__NodeElement__Group_2__2__Impl
            {
            pushFollow(FOLLOW_2);
            rule__NodeElement__Group_2__2__Impl();

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
    // $ANTLR end "rule__NodeElement__Group_2__2"


    // $ANTLR start "rule__NodeElement__Group_2__2__Impl"
    // InternalElkGraphJson.g:2035:1: rule__NodeElement__Group_2__2__Impl : ( ruleElkNodeChildren ) ;
    public final void rule__NodeElement__Group_2__2__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:2039:1: ( ( ruleElkNodeChildren ) )
            // InternalElkGraphJson.g:2040:1: ( ruleElkNodeChildren )
            {
            // InternalElkGraphJson.g:2040:1: ( ruleElkNodeChildren )
            // InternalElkGraphJson.g:2041:2: ruleElkNodeChildren
            {
             before(grammarAccess.getNodeElementAccess().getElkNodeChildrenParserRuleCall_2_2()); 
            pushFollow(FOLLOW_2);
            ruleElkNodeChildren();

            state._fsp--;

             after(grammarAccess.getNodeElementAccess().getElkNodeChildrenParserRuleCall_2_2()); 

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
    // $ANTLR end "rule__NodeElement__Group_2__2__Impl"


    // $ANTLR start "rule__NodeElement__Group_3__0"
    // InternalElkGraphJson.g:2051:1: rule__NodeElement__Group_3__0 : rule__NodeElement__Group_3__0__Impl rule__NodeElement__Group_3__1 ;
    public final void rule__NodeElement__Group_3__0() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:2055:1: ( rule__NodeElement__Group_3__0__Impl rule__NodeElement__Group_3__1 )
            // InternalElkGraphJson.g:2056:2: rule__NodeElement__Group_3__0__Impl rule__NodeElement__Group_3__1
            {
            pushFollow(FOLLOW_8);
            rule__NodeElement__Group_3__0__Impl();

            state._fsp--;

            pushFollow(FOLLOW_2);
            rule__NodeElement__Group_3__1();

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
    // $ANTLR end "rule__NodeElement__Group_3__0"


    // $ANTLR start "rule__NodeElement__Group_3__0__Impl"
    // InternalElkGraphJson.g:2063:1: rule__NodeElement__Group_3__0__Impl : ( ruleKeyPorts ) ;
    public final void rule__NodeElement__Group_3__0__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:2067:1: ( ( ruleKeyPorts ) )
            // InternalElkGraphJson.g:2068:1: ( ruleKeyPorts )
            {
            // InternalElkGraphJson.g:2068:1: ( ruleKeyPorts )
            // InternalElkGraphJson.g:2069:2: ruleKeyPorts
            {
             before(grammarAccess.getNodeElementAccess().getKeyPortsParserRuleCall_3_0()); 
            pushFollow(FOLLOW_2);
            ruleKeyPorts();

            state._fsp--;

             after(grammarAccess.getNodeElementAccess().getKeyPortsParserRuleCall_3_0()); 

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
    // $ANTLR end "rule__NodeElement__Group_3__0__Impl"


    // $ANTLR start "rule__NodeElement__Group_3__1"
    // InternalElkGraphJson.g:2078:1: rule__NodeElement__Group_3__1 : rule__NodeElement__Group_3__1__Impl rule__NodeElement__Group_3__2 ;
    public final void rule__NodeElement__Group_3__1() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:2082:1: ( rule__NodeElement__Group_3__1__Impl rule__NodeElement__Group_3__2 )
            // InternalElkGraphJson.g:2083:2: rule__NodeElement__Group_3__1__Impl rule__NodeElement__Group_3__2
            {
            pushFollow(FOLLOW_9);
            rule__NodeElement__Group_3__1__Impl();

            state._fsp--;

            pushFollow(FOLLOW_2);
            rule__NodeElement__Group_3__2();

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
    // $ANTLR end "rule__NodeElement__Group_3__1"


    // $ANTLR start "rule__NodeElement__Group_3__1__Impl"
    // InternalElkGraphJson.g:2090:1: rule__NodeElement__Group_3__1__Impl : ( ':' ) ;
    public final void rule__NodeElement__Group_3__1__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:2094:1: ( ( ':' ) )
            // InternalElkGraphJson.g:2095:1: ( ':' )
            {
            // InternalElkGraphJson.g:2095:1: ( ':' )
            // InternalElkGraphJson.g:2096:2: ':'
            {
             before(grammarAccess.getNodeElementAccess().getColonKeyword_3_1()); 
            match(input,61,FOLLOW_2); 
             after(grammarAccess.getNodeElementAccess().getColonKeyword_3_1()); 

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
    // $ANTLR end "rule__NodeElement__Group_3__1__Impl"


    // $ANTLR start "rule__NodeElement__Group_3__2"
    // InternalElkGraphJson.g:2105:1: rule__NodeElement__Group_3__2 : rule__NodeElement__Group_3__2__Impl ;
    public final void rule__NodeElement__Group_3__2() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:2109:1: ( rule__NodeElement__Group_3__2__Impl )
            // InternalElkGraphJson.g:2110:2: rule__NodeElement__Group_3__2__Impl
            {
            pushFollow(FOLLOW_2);
            rule__NodeElement__Group_3__2__Impl();

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
    // $ANTLR end "rule__NodeElement__Group_3__2"


    // $ANTLR start "rule__NodeElement__Group_3__2__Impl"
    // InternalElkGraphJson.g:2116:1: rule__NodeElement__Group_3__2__Impl : ( ruleElkNodePorts ) ;
    public final void rule__NodeElement__Group_3__2__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:2120:1: ( ( ruleElkNodePorts ) )
            // InternalElkGraphJson.g:2121:1: ( ruleElkNodePorts )
            {
            // InternalElkGraphJson.g:2121:1: ( ruleElkNodePorts )
            // InternalElkGraphJson.g:2122:2: ruleElkNodePorts
            {
             before(grammarAccess.getNodeElementAccess().getElkNodePortsParserRuleCall_3_2()); 
            pushFollow(FOLLOW_2);
            ruleElkNodePorts();

            state._fsp--;

             after(grammarAccess.getNodeElementAccess().getElkNodePortsParserRuleCall_3_2()); 

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
    // $ANTLR end "rule__NodeElement__Group_3__2__Impl"


    // $ANTLR start "rule__NodeElement__Group_4__0"
    // InternalElkGraphJson.g:2132:1: rule__NodeElement__Group_4__0 : rule__NodeElement__Group_4__0__Impl rule__NodeElement__Group_4__1 ;
    public final void rule__NodeElement__Group_4__0() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:2136:1: ( rule__NodeElement__Group_4__0__Impl rule__NodeElement__Group_4__1 )
            // InternalElkGraphJson.g:2137:2: rule__NodeElement__Group_4__0__Impl rule__NodeElement__Group_4__1
            {
            pushFollow(FOLLOW_8);
            rule__NodeElement__Group_4__0__Impl();

            state._fsp--;

            pushFollow(FOLLOW_2);
            rule__NodeElement__Group_4__1();

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
    // $ANTLR end "rule__NodeElement__Group_4__0"


    // $ANTLR start "rule__NodeElement__Group_4__0__Impl"
    // InternalElkGraphJson.g:2144:1: rule__NodeElement__Group_4__0__Impl : ( ruleKeyLabels ) ;
    public final void rule__NodeElement__Group_4__0__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:2148:1: ( ( ruleKeyLabels ) )
            // InternalElkGraphJson.g:2149:1: ( ruleKeyLabels )
            {
            // InternalElkGraphJson.g:2149:1: ( ruleKeyLabels )
            // InternalElkGraphJson.g:2150:2: ruleKeyLabels
            {
             before(grammarAccess.getNodeElementAccess().getKeyLabelsParserRuleCall_4_0()); 
            pushFollow(FOLLOW_2);
            ruleKeyLabels();

            state._fsp--;

             after(grammarAccess.getNodeElementAccess().getKeyLabelsParserRuleCall_4_0()); 

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
    // $ANTLR end "rule__NodeElement__Group_4__0__Impl"


    // $ANTLR start "rule__NodeElement__Group_4__1"
    // InternalElkGraphJson.g:2159:1: rule__NodeElement__Group_4__1 : rule__NodeElement__Group_4__1__Impl rule__NodeElement__Group_4__2 ;
    public final void rule__NodeElement__Group_4__1() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:2163:1: ( rule__NodeElement__Group_4__1__Impl rule__NodeElement__Group_4__2 )
            // InternalElkGraphJson.g:2164:2: rule__NodeElement__Group_4__1__Impl rule__NodeElement__Group_4__2
            {
            pushFollow(FOLLOW_9);
            rule__NodeElement__Group_4__1__Impl();

            state._fsp--;

            pushFollow(FOLLOW_2);
            rule__NodeElement__Group_4__2();

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
    // $ANTLR end "rule__NodeElement__Group_4__1"


    // $ANTLR start "rule__NodeElement__Group_4__1__Impl"
    // InternalElkGraphJson.g:2171:1: rule__NodeElement__Group_4__1__Impl : ( ':' ) ;
    public final void rule__NodeElement__Group_4__1__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:2175:1: ( ( ':' ) )
            // InternalElkGraphJson.g:2176:1: ( ':' )
            {
            // InternalElkGraphJson.g:2176:1: ( ':' )
            // InternalElkGraphJson.g:2177:2: ':'
            {
             before(grammarAccess.getNodeElementAccess().getColonKeyword_4_1()); 
            match(input,61,FOLLOW_2); 
             after(grammarAccess.getNodeElementAccess().getColonKeyword_4_1()); 

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
    // $ANTLR end "rule__NodeElement__Group_4__1__Impl"


    // $ANTLR start "rule__NodeElement__Group_4__2"
    // InternalElkGraphJson.g:2186:1: rule__NodeElement__Group_4__2 : rule__NodeElement__Group_4__2__Impl ;
    public final void rule__NodeElement__Group_4__2() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:2190:1: ( rule__NodeElement__Group_4__2__Impl )
            // InternalElkGraphJson.g:2191:2: rule__NodeElement__Group_4__2__Impl
            {
            pushFollow(FOLLOW_2);
            rule__NodeElement__Group_4__2__Impl();

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
    // $ANTLR end "rule__NodeElement__Group_4__2"


    // $ANTLR start "rule__NodeElement__Group_4__2__Impl"
    // InternalElkGraphJson.g:2197:1: rule__NodeElement__Group_4__2__Impl : ( ruleElkGraphElementLabels ) ;
    public final void rule__NodeElement__Group_4__2__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:2201:1: ( ( ruleElkGraphElementLabels ) )
            // InternalElkGraphJson.g:2202:1: ( ruleElkGraphElementLabels )
            {
            // InternalElkGraphJson.g:2202:1: ( ruleElkGraphElementLabels )
            // InternalElkGraphJson.g:2203:2: ruleElkGraphElementLabels
            {
             before(grammarAccess.getNodeElementAccess().getElkGraphElementLabelsParserRuleCall_4_2()); 
            pushFollow(FOLLOW_2);
            ruleElkGraphElementLabels();

            state._fsp--;

             after(grammarAccess.getNodeElementAccess().getElkGraphElementLabelsParserRuleCall_4_2()); 

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
    // $ANTLR end "rule__NodeElement__Group_4__2__Impl"


    // $ANTLR start "rule__NodeElement__Group_5__0"
    // InternalElkGraphJson.g:2213:1: rule__NodeElement__Group_5__0 : rule__NodeElement__Group_5__0__Impl rule__NodeElement__Group_5__1 ;
    public final void rule__NodeElement__Group_5__0() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:2217:1: ( rule__NodeElement__Group_5__0__Impl rule__NodeElement__Group_5__1 )
            // InternalElkGraphJson.g:2218:2: rule__NodeElement__Group_5__0__Impl rule__NodeElement__Group_5__1
            {
            pushFollow(FOLLOW_8);
            rule__NodeElement__Group_5__0__Impl();

            state._fsp--;

            pushFollow(FOLLOW_2);
            rule__NodeElement__Group_5__1();

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
    // $ANTLR end "rule__NodeElement__Group_5__0"


    // $ANTLR start "rule__NodeElement__Group_5__0__Impl"
    // InternalElkGraphJson.g:2225:1: rule__NodeElement__Group_5__0__Impl : ( ruleKeyEdges ) ;
    public final void rule__NodeElement__Group_5__0__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:2229:1: ( ( ruleKeyEdges ) )
            // InternalElkGraphJson.g:2230:1: ( ruleKeyEdges )
            {
            // InternalElkGraphJson.g:2230:1: ( ruleKeyEdges )
            // InternalElkGraphJson.g:2231:2: ruleKeyEdges
            {
             before(grammarAccess.getNodeElementAccess().getKeyEdgesParserRuleCall_5_0()); 
            pushFollow(FOLLOW_2);
            ruleKeyEdges();

            state._fsp--;

             after(grammarAccess.getNodeElementAccess().getKeyEdgesParserRuleCall_5_0()); 

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
    // $ANTLR end "rule__NodeElement__Group_5__0__Impl"


    // $ANTLR start "rule__NodeElement__Group_5__1"
    // InternalElkGraphJson.g:2240:1: rule__NodeElement__Group_5__1 : rule__NodeElement__Group_5__1__Impl rule__NodeElement__Group_5__2 ;
    public final void rule__NodeElement__Group_5__1() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:2244:1: ( rule__NodeElement__Group_5__1__Impl rule__NodeElement__Group_5__2 )
            // InternalElkGraphJson.g:2245:2: rule__NodeElement__Group_5__1__Impl rule__NodeElement__Group_5__2
            {
            pushFollow(FOLLOW_9);
            rule__NodeElement__Group_5__1__Impl();

            state._fsp--;

            pushFollow(FOLLOW_2);
            rule__NodeElement__Group_5__2();

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
    // $ANTLR end "rule__NodeElement__Group_5__1"


    // $ANTLR start "rule__NodeElement__Group_5__1__Impl"
    // InternalElkGraphJson.g:2252:1: rule__NodeElement__Group_5__1__Impl : ( ':' ) ;
    public final void rule__NodeElement__Group_5__1__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:2256:1: ( ( ':' ) )
            // InternalElkGraphJson.g:2257:1: ( ':' )
            {
            // InternalElkGraphJson.g:2257:1: ( ':' )
            // InternalElkGraphJson.g:2258:2: ':'
            {
             before(grammarAccess.getNodeElementAccess().getColonKeyword_5_1()); 
            match(input,61,FOLLOW_2); 
             after(grammarAccess.getNodeElementAccess().getColonKeyword_5_1()); 

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
    // $ANTLR end "rule__NodeElement__Group_5__1__Impl"


    // $ANTLR start "rule__NodeElement__Group_5__2"
    // InternalElkGraphJson.g:2267:1: rule__NodeElement__Group_5__2 : rule__NodeElement__Group_5__2__Impl ;
    public final void rule__NodeElement__Group_5__2() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:2271:1: ( rule__NodeElement__Group_5__2__Impl )
            // InternalElkGraphJson.g:2272:2: rule__NodeElement__Group_5__2__Impl
            {
            pushFollow(FOLLOW_2);
            rule__NodeElement__Group_5__2__Impl();

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
    // $ANTLR end "rule__NodeElement__Group_5__2"


    // $ANTLR start "rule__NodeElement__Group_5__2__Impl"
    // InternalElkGraphJson.g:2278:1: rule__NodeElement__Group_5__2__Impl : ( ruleElkNodeEdges ) ;
    public final void rule__NodeElement__Group_5__2__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:2282:1: ( ( ruleElkNodeEdges ) )
            // InternalElkGraphJson.g:2283:1: ( ruleElkNodeEdges )
            {
            // InternalElkGraphJson.g:2283:1: ( ruleElkNodeEdges )
            // InternalElkGraphJson.g:2284:2: ruleElkNodeEdges
            {
             before(grammarAccess.getNodeElementAccess().getElkNodeEdgesParserRuleCall_5_2()); 
            pushFollow(FOLLOW_2);
            ruleElkNodeEdges();

            state._fsp--;

             after(grammarAccess.getNodeElementAccess().getElkNodeEdgesParserRuleCall_5_2()); 

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
    // $ANTLR end "rule__NodeElement__Group_5__2__Impl"


    // $ANTLR start "rule__NodeElement__Group_6__0"
    // InternalElkGraphJson.g:2294:1: rule__NodeElement__Group_6__0 : rule__NodeElement__Group_6__0__Impl rule__NodeElement__Group_6__1 ;
    public final void rule__NodeElement__Group_6__0() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:2298:1: ( rule__NodeElement__Group_6__0__Impl rule__NodeElement__Group_6__1 )
            // InternalElkGraphJson.g:2299:2: rule__NodeElement__Group_6__0__Impl rule__NodeElement__Group_6__1
            {
            pushFollow(FOLLOW_8);
            rule__NodeElement__Group_6__0__Impl();

            state._fsp--;

            pushFollow(FOLLOW_2);
            rule__NodeElement__Group_6__1();

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
    // $ANTLR end "rule__NodeElement__Group_6__0"


    // $ANTLR start "rule__NodeElement__Group_6__0__Impl"
    // InternalElkGraphJson.g:2306:1: rule__NodeElement__Group_6__0__Impl : ( ruleKeyLayoutOptions ) ;
    public final void rule__NodeElement__Group_6__0__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:2310:1: ( ( ruleKeyLayoutOptions ) )
            // InternalElkGraphJson.g:2311:1: ( ruleKeyLayoutOptions )
            {
            // InternalElkGraphJson.g:2311:1: ( ruleKeyLayoutOptions )
            // InternalElkGraphJson.g:2312:2: ruleKeyLayoutOptions
            {
             before(grammarAccess.getNodeElementAccess().getKeyLayoutOptionsParserRuleCall_6_0()); 
            pushFollow(FOLLOW_2);
            ruleKeyLayoutOptions();

            state._fsp--;

             after(grammarAccess.getNodeElementAccess().getKeyLayoutOptionsParserRuleCall_6_0()); 

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
    // $ANTLR end "rule__NodeElement__Group_6__0__Impl"


    // $ANTLR start "rule__NodeElement__Group_6__1"
    // InternalElkGraphJson.g:2321:1: rule__NodeElement__Group_6__1 : rule__NodeElement__Group_6__1__Impl rule__NodeElement__Group_6__2 ;
    public final void rule__NodeElement__Group_6__1() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:2325:1: ( rule__NodeElement__Group_6__1__Impl rule__NodeElement__Group_6__2 )
            // InternalElkGraphJson.g:2326:2: rule__NodeElement__Group_6__1__Impl rule__NodeElement__Group_6__2
            {
            pushFollow(FOLLOW_3);
            rule__NodeElement__Group_6__1__Impl();

            state._fsp--;

            pushFollow(FOLLOW_2);
            rule__NodeElement__Group_6__2();

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
    // $ANTLR end "rule__NodeElement__Group_6__1"


    // $ANTLR start "rule__NodeElement__Group_6__1__Impl"
    // InternalElkGraphJson.g:2333:1: rule__NodeElement__Group_6__1__Impl : ( ':' ) ;
    public final void rule__NodeElement__Group_6__1__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:2337:1: ( ( ':' ) )
            // InternalElkGraphJson.g:2338:1: ( ':' )
            {
            // InternalElkGraphJson.g:2338:1: ( ':' )
            // InternalElkGraphJson.g:2339:2: ':'
            {
             before(grammarAccess.getNodeElementAccess().getColonKeyword_6_1()); 
            match(input,61,FOLLOW_2); 
             after(grammarAccess.getNodeElementAccess().getColonKeyword_6_1()); 

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
    // $ANTLR end "rule__NodeElement__Group_6__1__Impl"


    // $ANTLR start "rule__NodeElement__Group_6__2"
    // InternalElkGraphJson.g:2348:1: rule__NodeElement__Group_6__2 : rule__NodeElement__Group_6__2__Impl ;
    public final void rule__NodeElement__Group_6__2() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:2352:1: ( rule__NodeElement__Group_6__2__Impl )
            // InternalElkGraphJson.g:2353:2: rule__NodeElement__Group_6__2__Impl
            {
            pushFollow(FOLLOW_2);
            rule__NodeElement__Group_6__2__Impl();

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
    // $ANTLR end "rule__NodeElement__Group_6__2"


    // $ANTLR start "rule__NodeElement__Group_6__2__Impl"
    // InternalElkGraphJson.g:2359:1: rule__NodeElement__Group_6__2__Impl : ( ruleElkGraphElementProperties ) ;
    public final void rule__NodeElement__Group_6__2__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:2363:1: ( ( ruleElkGraphElementProperties ) )
            // InternalElkGraphJson.g:2364:1: ( ruleElkGraphElementProperties )
            {
            // InternalElkGraphJson.g:2364:1: ( ruleElkGraphElementProperties )
            // InternalElkGraphJson.g:2365:2: ruleElkGraphElementProperties
            {
             before(grammarAccess.getNodeElementAccess().getElkGraphElementPropertiesParserRuleCall_6_2()); 
            pushFollow(FOLLOW_2);
            ruleElkGraphElementProperties();

            state._fsp--;

             after(grammarAccess.getNodeElementAccess().getElkGraphElementPropertiesParserRuleCall_6_2()); 

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
    // $ANTLR end "rule__NodeElement__Group_6__2__Impl"


    // $ANTLR start "rule__ElkPort__Group__0"
    // InternalElkGraphJson.g:2375:1: rule__ElkPort__Group__0 : rule__ElkPort__Group__0__Impl rule__ElkPort__Group__1 ;
    public final void rule__ElkPort__Group__0() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:2379:1: ( rule__ElkPort__Group__0__Impl rule__ElkPort__Group__1 )
            // InternalElkGraphJson.g:2380:2: rule__ElkPort__Group__0__Impl rule__ElkPort__Group__1
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
    // InternalElkGraphJson.g:2387:1: rule__ElkPort__Group__0__Impl : ( '{' ) ;
    public final void rule__ElkPort__Group__0__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:2391:1: ( ( '{' ) )
            // InternalElkGraphJson.g:2392:1: ( '{' )
            {
            // InternalElkGraphJson.g:2392:1: ( '{' )
            // InternalElkGraphJson.g:2393:2: '{'
            {
             before(grammarAccess.getElkPortAccess().getLeftCurlyBracketKeyword_0()); 
            match(input,58,FOLLOW_2); 
             after(grammarAccess.getElkPortAccess().getLeftCurlyBracketKeyword_0()); 

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
    // InternalElkGraphJson.g:2402:1: rule__ElkPort__Group__1 : rule__ElkPort__Group__1__Impl rule__ElkPort__Group__2 ;
    public final void rule__ElkPort__Group__1() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:2406:1: ( rule__ElkPort__Group__1__Impl rule__ElkPort__Group__2 )
            // InternalElkGraphJson.g:2407:2: rule__ElkPort__Group__1__Impl rule__ElkPort__Group__2
            {
            pushFollow(FOLLOW_10);
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
    // InternalElkGraphJson.g:2414:1: rule__ElkPort__Group__1__Impl : ( rulePortElement ) ;
    public final void rule__ElkPort__Group__1__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:2418:1: ( ( rulePortElement ) )
            // InternalElkGraphJson.g:2419:1: ( rulePortElement )
            {
            // InternalElkGraphJson.g:2419:1: ( rulePortElement )
            // InternalElkGraphJson.g:2420:2: rulePortElement
            {
             before(grammarAccess.getElkPortAccess().getPortElementParserRuleCall_1()); 
            pushFollow(FOLLOW_2);
            rulePortElement();

            state._fsp--;

             after(grammarAccess.getElkPortAccess().getPortElementParserRuleCall_1()); 

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
    // InternalElkGraphJson.g:2429:1: rule__ElkPort__Group__2 : rule__ElkPort__Group__2__Impl rule__ElkPort__Group__3 ;
    public final void rule__ElkPort__Group__2() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:2433:1: ( rule__ElkPort__Group__2__Impl rule__ElkPort__Group__3 )
            // InternalElkGraphJson.g:2434:2: rule__ElkPort__Group__2__Impl rule__ElkPort__Group__3
            {
            pushFollow(FOLLOW_10);
            rule__ElkPort__Group__2__Impl();

            state._fsp--;

            pushFollow(FOLLOW_2);
            rule__ElkPort__Group__3();

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
    // InternalElkGraphJson.g:2441:1: rule__ElkPort__Group__2__Impl : ( ( rule__ElkPort__Group_2__0 )* ) ;
    public final void rule__ElkPort__Group__2__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:2445:1: ( ( ( rule__ElkPort__Group_2__0 )* ) )
            // InternalElkGraphJson.g:2446:1: ( ( rule__ElkPort__Group_2__0 )* )
            {
            // InternalElkGraphJson.g:2446:1: ( ( rule__ElkPort__Group_2__0 )* )
            // InternalElkGraphJson.g:2447:2: ( rule__ElkPort__Group_2__0 )*
            {
             before(grammarAccess.getElkPortAccess().getGroup_2()); 
            // InternalElkGraphJson.g:2448:2: ( rule__ElkPort__Group_2__0 )*
            loop29:
            do {
                int alt29=2;
                int LA29_0 = input.LA(1);

                if ( (LA29_0==59) ) {
                    int LA29_1 = input.LA(2);

                    if ( ((LA29_1>=RULE_STRING && LA29_1<=RULE_ID)||(LA29_1>=22 && LA29_1<=24)||(LA29_1>=28 && LA29_1<=48)) ) {
                        alt29=1;
                    }


                }


                switch (alt29) {
            	case 1 :
            	    // InternalElkGraphJson.g:2448:3: rule__ElkPort__Group_2__0
            	    {
            	    pushFollow(FOLLOW_6);
            	    rule__ElkPort__Group_2__0();

            	    state._fsp--;


            	    }
            	    break;

            	default :
            	    break loop29;
                }
            } while (true);

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


    // $ANTLR start "rule__ElkPort__Group__3"
    // InternalElkGraphJson.g:2456:1: rule__ElkPort__Group__3 : rule__ElkPort__Group__3__Impl rule__ElkPort__Group__4 ;
    public final void rule__ElkPort__Group__3() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:2460:1: ( rule__ElkPort__Group__3__Impl rule__ElkPort__Group__4 )
            // InternalElkGraphJson.g:2461:2: rule__ElkPort__Group__3__Impl rule__ElkPort__Group__4
            {
            pushFollow(FOLLOW_10);
            rule__ElkPort__Group__3__Impl();

            state._fsp--;

            pushFollow(FOLLOW_2);
            rule__ElkPort__Group__4();

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
    // $ANTLR end "rule__ElkPort__Group__3"


    // $ANTLR start "rule__ElkPort__Group__3__Impl"
    // InternalElkGraphJson.g:2468:1: rule__ElkPort__Group__3__Impl : ( ( ',' )? ) ;
    public final void rule__ElkPort__Group__3__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:2472:1: ( ( ( ',' )? ) )
            // InternalElkGraphJson.g:2473:1: ( ( ',' )? )
            {
            // InternalElkGraphJson.g:2473:1: ( ( ',' )? )
            // InternalElkGraphJson.g:2474:2: ( ',' )?
            {
             before(grammarAccess.getElkPortAccess().getCommaKeyword_3()); 
            // InternalElkGraphJson.g:2475:2: ( ',' )?
            int alt30=2;
            int LA30_0 = input.LA(1);

            if ( (LA30_0==59) ) {
                alt30=1;
            }
            switch (alt30) {
                case 1 :
                    // InternalElkGraphJson.g:2475:3: ','
                    {
                    match(input,59,FOLLOW_2); 

                    }
                    break;

            }

             after(grammarAccess.getElkPortAccess().getCommaKeyword_3()); 

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
    // $ANTLR end "rule__ElkPort__Group__3__Impl"


    // $ANTLR start "rule__ElkPort__Group__4"
    // InternalElkGraphJson.g:2483:1: rule__ElkPort__Group__4 : rule__ElkPort__Group__4__Impl ;
    public final void rule__ElkPort__Group__4() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:2487:1: ( rule__ElkPort__Group__4__Impl )
            // InternalElkGraphJson.g:2488:2: rule__ElkPort__Group__4__Impl
            {
            pushFollow(FOLLOW_2);
            rule__ElkPort__Group__4__Impl();

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
    // $ANTLR end "rule__ElkPort__Group__4"


    // $ANTLR start "rule__ElkPort__Group__4__Impl"
    // InternalElkGraphJson.g:2494:1: rule__ElkPort__Group__4__Impl : ( '}' ) ;
    public final void rule__ElkPort__Group__4__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:2498:1: ( ( '}' ) )
            // InternalElkGraphJson.g:2499:1: ( '}' )
            {
            // InternalElkGraphJson.g:2499:1: ( '}' )
            // InternalElkGraphJson.g:2500:2: '}'
            {
             before(grammarAccess.getElkPortAccess().getRightCurlyBracketKeyword_4()); 
            match(input,60,FOLLOW_2); 
             after(grammarAccess.getElkPortAccess().getRightCurlyBracketKeyword_4()); 

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
    // $ANTLR end "rule__ElkPort__Group__4__Impl"


    // $ANTLR start "rule__ElkPort__Group_2__0"
    // InternalElkGraphJson.g:2510:1: rule__ElkPort__Group_2__0 : rule__ElkPort__Group_2__0__Impl rule__ElkPort__Group_2__1 ;
    public final void rule__ElkPort__Group_2__0() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:2514:1: ( rule__ElkPort__Group_2__0__Impl rule__ElkPort__Group_2__1 )
            // InternalElkGraphJson.g:2515:2: rule__ElkPort__Group_2__0__Impl rule__ElkPort__Group_2__1
            {
            pushFollow(FOLLOW_7);
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
    // InternalElkGraphJson.g:2522:1: rule__ElkPort__Group_2__0__Impl : ( ',' ) ;
    public final void rule__ElkPort__Group_2__0__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:2526:1: ( ( ',' ) )
            // InternalElkGraphJson.g:2527:1: ( ',' )
            {
            // InternalElkGraphJson.g:2527:1: ( ',' )
            // InternalElkGraphJson.g:2528:2: ','
            {
             before(grammarAccess.getElkPortAccess().getCommaKeyword_2_0()); 
            match(input,59,FOLLOW_2); 
             after(grammarAccess.getElkPortAccess().getCommaKeyword_2_0()); 

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
    // InternalElkGraphJson.g:2537:1: rule__ElkPort__Group_2__1 : rule__ElkPort__Group_2__1__Impl ;
    public final void rule__ElkPort__Group_2__1() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:2541:1: ( rule__ElkPort__Group_2__1__Impl )
            // InternalElkGraphJson.g:2542:2: rule__ElkPort__Group_2__1__Impl
            {
            pushFollow(FOLLOW_2);
            rule__ElkPort__Group_2__1__Impl();

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
    // InternalElkGraphJson.g:2548:1: rule__ElkPort__Group_2__1__Impl : ( rulePortElement ) ;
    public final void rule__ElkPort__Group_2__1__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:2552:1: ( ( rulePortElement ) )
            // InternalElkGraphJson.g:2553:1: ( rulePortElement )
            {
            // InternalElkGraphJson.g:2553:1: ( rulePortElement )
            // InternalElkGraphJson.g:2554:2: rulePortElement
            {
             before(grammarAccess.getElkPortAccess().getPortElementParserRuleCall_2_1()); 
            pushFollow(FOLLOW_2);
            rulePortElement();

            state._fsp--;

             after(grammarAccess.getElkPortAccess().getPortElementParserRuleCall_2_1()); 

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


    // $ANTLR start "rule__PortElement__Group_2__0"
    // InternalElkGraphJson.g:2564:1: rule__PortElement__Group_2__0 : rule__PortElement__Group_2__0__Impl rule__PortElement__Group_2__1 ;
    public final void rule__PortElement__Group_2__0() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:2568:1: ( rule__PortElement__Group_2__0__Impl rule__PortElement__Group_2__1 )
            // InternalElkGraphJson.g:2569:2: rule__PortElement__Group_2__0__Impl rule__PortElement__Group_2__1
            {
            pushFollow(FOLLOW_8);
            rule__PortElement__Group_2__0__Impl();

            state._fsp--;

            pushFollow(FOLLOW_2);
            rule__PortElement__Group_2__1();

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
    // $ANTLR end "rule__PortElement__Group_2__0"


    // $ANTLR start "rule__PortElement__Group_2__0__Impl"
    // InternalElkGraphJson.g:2576:1: rule__PortElement__Group_2__0__Impl : ( ruleKeyLabels ) ;
    public final void rule__PortElement__Group_2__0__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:2580:1: ( ( ruleKeyLabels ) )
            // InternalElkGraphJson.g:2581:1: ( ruleKeyLabels )
            {
            // InternalElkGraphJson.g:2581:1: ( ruleKeyLabels )
            // InternalElkGraphJson.g:2582:2: ruleKeyLabels
            {
             before(grammarAccess.getPortElementAccess().getKeyLabelsParserRuleCall_2_0()); 
            pushFollow(FOLLOW_2);
            ruleKeyLabels();

            state._fsp--;

             after(grammarAccess.getPortElementAccess().getKeyLabelsParserRuleCall_2_0()); 

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
    // $ANTLR end "rule__PortElement__Group_2__0__Impl"


    // $ANTLR start "rule__PortElement__Group_2__1"
    // InternalElkGraphJson.g:2591:1: rule__PortElement__Group_2__1 : rule__PortElement__Group_2__1__Impl rule__PortElement__Group_2__2 ;
    public final void rule__PortElement__Group_2__1() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:2595:1: ( rule__PortElement__Group_2__1__Impl rule__PortElement__Group_2__2 )
            // InternalElkGraphJson.g:2596:2: rule__PortElement__Group_2__1__Impl rule__PortElement__Group_2__2
            {
            pushFollow(FOLLOW_9);
            rule__PortElement__Group_2__1__Impl();

            state._fsp--;

            pushFollow(FOLLOW_2);
            rule__PortElement__Group_2__2();

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
    // $ANTLR end "rule__PortElement__Group_2__1"


    // $ANTLR start "rule__PortElement__Group_2__1__Impl"
    // InternalElkGraphJson.g:2603:1: rule__PortElement__Group_2__1__Impl : ( ':' ) ;
    public final void rule__PortElement__Group_2__1__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:2607:1: ( ( ':' ) )
            // InternalElkGraphJson.g:2608:1: ( ':' )
            {
            // InternalElkGraphJson.g:2608:1: ( ':' )
            // InternalElkGraphJson.g:2609:2: ':'
            {
             before(grammarAccess.getPortElementAccess().getColonKeyword_2_1()); 
            match(input,61,FOLLOW_2); 
             after(grammarAccess.getPortElementAccess().getColonKeyword_2_1()); 

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
    // $ANTLR end "rule__PortElement__Group_2__1__Impl"


    // $ANTLR start "rule__PortElement__Group_2__2"
    // InternalElkGraphJson.g:2618:1: rule__PortElement__Group_2__2 : rule__PortElement__Group_2__2__Impl ;
    public final void rule__PortElement__Group_2__2() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:2622:1: ( rule__PortElement__Group_2__2__Impl )
            // InternalElkGraphJson.g:2623:2: rule__PortElement__Group_2__2__Impl
            {
            pushFollow(FOLLOW_2);
            rule__PortElement__Group_2__2__Impl();

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
    // $ANTLR end "rule__PortElement__Group_2__2"


    // $ANTLR start "rule__PortElement__Group_2__2__Impl"
    // InternalElkGraphJson.g:2629:1: rule__PortElement__Group_2__2__Impl : ( ruleElkGraphElementLabels ) ;
    public final void rule__PortElement__Group_2__2__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:2633:1: ( ( ruleElkGraphElementLabels ) )
            // InternalElkGraphJson.g:2634:1: ( ruleElkGraphElementLabels )
            {
            // InternalElkGraphJson.g:2634:1: ( ruleElkGraphElementLabels )
            // InternalElkGraphJson.g:2635:2: ruleElkGraphElementLabels
            {
             before(grammarAccess.getPortElementAccess().getElkGraphElementLabelsParserRuleCall_2_2()); 
            pushFollow(FOLLOW_2);
            ruleElkGraphElementLabels();

            state._fsp--;

             after(grammarAccess.getPortElementAccess().getElkGraphElementLabelsParserRuleCall_2_2()); 

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
    // $ANTLR end "rule__PortElement__Group_2__2__Impl"


    // $ANTLR start "rule__PortElement__Group_3__0"
    // InternalElkGraphJson.g:2645:1: rule__PortElement__Group_3__0 : rule__PortElement__Group_3__0__Impl rule__PortElement__Group_3__1 ;
    public final void rule__PortElement__Group_3__0() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:2649:1: ( rule__PortElement__Group_3__0__Impl rule__PortElement__Group_3__1 )
            // InternalElkGraphJson.g:2650:2: rule__PortElement__Group_3__0__Impl rule__PortElement__Group_3__1
            {
            pushFollow(FOLLOW_8);
            rule__PortElement__Group_3__0__Impl();

            state._fsp--;

            pushFollow(FOLLOW_2);
            rule__PortElement__Group_3__1();

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
    // $ANTLR end "rule__PortElement__Group_3__0"


    // $ANTLR start "rule__PortElement__Group_3__0__Impl"
    // InternalElkGraphJson.g:2657:1: rule__PortElement__Group_3__0__Impl : ( ruleKeyLayoutOptions ) ;
    public final void rule__PortElement__Group_3__0__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:2661:1: ( ( ruleKeyLayoutOptions ) )
            // InternalElkGraphJson.g:2662:1: ( ruleKeyLayoutOptions )
            {
            // InternalElkGraphJson.g:2662:1: ( ruleKeyLayoutOptions )
            // InternalElkGraphJson.g:2663:2: ruleKeyLayoutOptions
            {
             before(grammarAccess.getPortElementAccess().getKeyLayoutOptionsParserRuleCall_3_0()); 
            pushFollow(FOLLOW_2);
            ruleKeyLayoutOptions();

            state._fsp--;

             after(grammarAccess.getPortElementAccess().getKeyLayoutOptionsParserRuleCall_3_0()); 

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
    // $ANTLR end "rule__PortElement__Group_3__0__Impl"


    // $ANTLR start "rule__PortElement__Group_3__1"
    // InternalElkGraphJson.g:2672:1: rule__PortElement__Group_3__1 : rule__PortElement__Group_3__1__Impl rule__PortElement__Group_3__2 ;
    public final void rule__PortElement__Group_3__1() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:2676:1: ( rule__PortElement__Group_3__1__Impl rule__PortElement__Group_3__2 )
            // InternalElkGraphJson.g:2677:2: rule__PortElement__Group_3__1__Impl rule__PortElement__Group_3__2
            {
            pushFollow(FOLLOW_3);
            rule__PortElement__Group_3__1__Impl();

            state._fsp--;

            pushFollow(FOLLOW_2);
            rule__PortElement__Group_3__2();

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
    // $ANTLR end "rule__PortElement__Group_3__1"


    // $ANTLR start "rule__PortElement__Group_3__1__Impl"
    // InternalElkGraphJson.g:2684:1: rule__PortElement__Group_3__1__Impl : ( ':' ) ;
    public final void rule__PortElement__Group_3__1__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:2688:1: ( ( ':' ) )
            // InternalElkGraphJson.g:2689:1: ( ':' )
            {
            // InternalElkGraphJson.g:2689:1: ( ':' )
            // InternalElkGraphJson.g:2690:2: ':'
            {
             before(grammarAccess.getPortElementAccess().getColonKeyword_3_1()); 
            match(input,61,FOLLOW_2); 
             after(grammarAccess.getPortElementAccess().getColonKeyword_3_1()); 

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
    // $ANTLR end "rule__PortElement__Group_3__1__Impl"


    // $ANTLR start "rule__PortElement__Group_3__2"
    // InternalElkGraphJson.g:2699:1: rule__PortElement__Group_3__2 : rule__PortElement__Group_3__2__Impl ;
    public final void rule__PortElement__Group_3__2() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:2703:1: ( rule__PortElement__Group_3__2__Impl )
            // InternalElkGraphJson.g:2704:2: rule__PortElement__Group_3__2__Impl
            {
            pushFollow(FOLLOW_2);
            rule__PortElement__Group_3__2__Impl();

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
    // $ANTLR end "rule__PortElement__Group_3__2"


    // $ANTLR start "rule__PortElement__Group_3__2__Impl"
    // InternalElkGraphJson.g:2710:1: rule__PortElement__Group_3__2__Impl : ( ruleElkGraphElementProperties ) ;
    public final void rule__PortElement__Group_3__2__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:2714:1: ( ( ruleElkGraphElementProperties ) )
            // InternalElkGraphJson.g:2715:1: ( ruleElkGraphElementProperties )
            {
            // InternalElkGraphJson.g:2715:1: ( ruleElkGraphElementProperties )
            // InternalElkGraphJson.g:2716:2: ruleElkGraphElementProperties
            {
             before(grammarAccess.getPortElementAccess().getElkGraphElementPropertiesParserRuleCall_3_2()); 
            pushFollow(FOLLOW_2);
            ruleElkGraphElementProperties();

            state._fsp--;

             after(grammarAccess.getPortElementAccess().getElkGraphElementPropertiesParserRuleCall_3_2()); 

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
    // $ANTLR end "rule__PortElement__Group_3__2__Impl"


    // $ANTLR start "rule__ElkLabel__Group__0"
    // InternalElkGraphJson.g:2726:1: rule__ElkLabel__Group__0 : rule__ElkLabel__Group__0__Impl rule__ElkLabel__Group__1 ;
    public final void rule__ElkLabel__Group__0() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:2730:1: ( rule__ElkLabel__Group__0__Impl rule__ElkLabel__Group__1 )
            // InternalElkGraphJson.g:2731:2: rule__ElkLabel__Group__0__Impl rule__ElkLabel__Group__1
            {
            pushFollow(FOLLOW_11);
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
    // InternalElkGraphJson.g:2738:1: rule__ElkLabel__Group__0__Impl : ( '{' ) ;
    public final void rule__ElkLabel__Group__0__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:2742:1: ( ( '{' ) )
            // InternalElkGraphJson.g:2743:1: ( '{' )
            {
            // InternalElkGraphJson.g:2743:1: ( '{' )
            // InternalElkGraphJson.g:2744:2: '{'
            {
             before(grammarAccess.getElkLabelAccess().getLeftCurlyBracketKeyword_0()); 
            match(input,58,FOLLOW_2); 
             after(grammarAccess.getElkLabelAccess().getLeftCurlyBracketKeyword_0()); 

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
    // InternalElkGraphJson.g:2753:1: rule__ElkLabel__Group__1 : rule__ElkLabel__Group__1__Impl rule__ElkLabel__Group__2 ;
    public final void rule__ElkLabel__Group__1() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:2757:1: ( rule__ElkLabel__Group__1__Impl rule__ElkLabel__Group__2 )
            // InternalElkGraphJson.g:2758:2: rule__ElkLabel__Group__1__Impl rule__ElkLabel__Group__2
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
    // InternalElkGraphJson.g:2765:1: rule__ElkLabel__Group__1__Impl : ( ruleLabelElement ) ;
    public final void rule__ElkLabel__Group__1__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:2769:1: ( ( ruleLabelElement ) )
            // InternalElkGraphJson.g:2770:1: ( ruleLabelElement )
            {
            // InternalElkGraphJson.g:2770:1: ( ruleLabelElement )
            // InternalElkGraphJson.g:2771:2: ruleLabelElement
            {
             before(grammarAccess.getElkLabelAccess().getLabelElementParserRuleCall_1()); 
            pushFollow(FOLLOW_2);
            ruleLabelElement();

            state._fsp--;

             after(grammarAccess.getElkLabelAccess().getLabelElementParserRuleCall_1()); 

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
    // InternalElkGraphJson.g:2780:1: rule__ElkLabel__Group__2 : rule__ElkLabel__Group__2__Impl rule__ElkLabel__Group__3 ;
    public final void rule__ElkLabel__Group__2() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:2784:1: ( rule__ElkLabel__Group__2__Impl rule__ElkLabel__Group__3 )
            // InternalElkGraphJson.g:2785:2: rule__ElkLabel__Group__2__Impl rule__ElkLabel__Group__3
            {
            pushFollow(FOLLOW_10);
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
    // InternalElkGraphJson.g:2792:1: rule__ElkLabel__Group__2__Impl : ( ( rule__ElkLabel__Group_2__0 )* ) ;
    public final void rule__ElkLabel__Group__2__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:2796:1: ( ( ( rule__ElkLabel__Group_2__0 )* ) )
            // InternalElkGraphJson.g:2797:1: ( ( rule__ElkLabel__Group_2__0 )* )
            {
            // InternalElkGraphJson.g:2797:1: ( ( rule__ElkLabel__Group_2__0 )* )
            // InternalElkGraphJson.g:2798:2: ( rule__ElkLabel__Group_2__0 )*
            {
             before(grammarAccess.getElkLabelAccess().getGroup_2()); 
            // InternalElkGraphJson.g:2799:2: ( rule__ElkLabel__Group_2__0 )*
            loop31:
            do {
                int alt31=2;
                int LA31_0 = input.LA(1);

                if ( (LA31_0==59) ) {
                    int LA31_1 = input.LA(2);

                    if ( ((LA31_1>=RULE_STRING && LA31_1<=RULE_ID)||(LA31_1>=22 && LA31_1<=24)||(LA31_1>=28 && LA31_1<=48)||(LA31_1>=55 && LA31_1<=57)) ) {
                        alt31=1;
                    }


                }


                switch (alt31) {
            	case 1 :
            	    // InternalElkGraphJson.g:2799:3: rule__ElkLabel__Group_2__0
            	    {
            	    pushFollow(FOLLOW_6);
            	    rule__ElkLabel__Group_2__0();

            	    state._fsp--;


            	    }
            	    break;

            	default :
            	    break loop31;
                }
            } while (true);

             after(grammarAccess.getElkLabelAccess().getGroup_2()); 

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
    // InternalElkGraphJson.g:2807:1: rule__ElkLabel__Group__3 : rule__ElkLabel__Group__3__Impl rule__ElkLabel__Group__4 ;
    public final void rule__ElkLabel__Group__3() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:2811:1: ( rule__ElkLabel__Group__3__Impl rule__ElkLabel__Group__4 )
            // InternalElkGraphJson.g:2812:2: rule__ElkLabel__Group__3__Impl rule__ElkLabel__Group__4
            {
            pushFollow(FOLLOW_10);
            rule__ElkLabel__Group__3__Impl();

            state._fsp--;

            pushFollow(FOLLOW_2);
            rule__ElkLabel__Group__4();

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
    // InternalElkGraphJson.g:2819:1: rule__ElkLabel__Group__3__Impl : ( ( ',' )? ) ;
    public final void rule__ElkLabel__Group__3__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:2823:1: ( ( ( ',' )? ) )
            // InternalElkGraphJson.g:2824:1: ( ( ',' )? )
            {
            // InternalElkGraphJson.g:2824:1: ( ( ',' )? )
            // InternalElkGraphJson.g:2825:2: ( ',' )?
            {
             before(grammarAccess.getElkLabelAccess().getCommaKeyword_3()); 
            // InternalElkGraphJson.g:2826:2: ( ',' )?
            int alt32=2;
            int LA32_0 = input.LA(1);

            if ( (LA32_0==59) ) {
                alt32=1;
            }
            switch (alt32) {
                case 1 :
                    // InternalElkGraphJson.g:2826:3: ','
                    {
                    match(input,59,FOLLOW_2); 

                    }
                    break;

            }

             after(grammarAccess.getElkLabelAccess().getCommaKeyword_3()); 

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


    // $ANTLR start "rule__ElkLabel__Group__4"
    // InternalElkGraphJson.g:2834:1: rule__ElkLabel__Group__4 : rule__ElkLabel__Group__4__Impl ;
    public final void rule__ElkLabel__Group__4() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:2838:1: ( rule__ElkLabel__Group__4__Impl )
            // InternalElkGraphJson.g:2839:2: rule__ElkLabel__Group__4__Impl
            {
            pushFollow(FOLLOW_2);
            rule__ElkLabel__Group__4__Impl();

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
    // $ANTLR end "rule__ElkLabel__Group__4"


    // $ANTLR start "rule__ElkLabel__Group__4__Impl"
    // InternalElkGraphJson.g:2845:1: rule__ElkLabel__Group__4__Impl : ( '}' ) ;
    public final void rule__ElkLabel__Group__4__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:2849:1: ( ( '}' ) )
            // InternalElkGraphJson.g:2850:1: ( '}' )
            {
            // InternalElkGraphJson.g:2850:1: ( '}' )
            // InternalElkGraphJson.g:2851:2: '}'
            {
             before(grammarAccess.getElkLabelAccess().getRightCurlyBracketKeyword_4()); 
            match(input,60,FOLLOW_2); 
             after(grammarAccess.getElkLabelAccess().getRightCurlyBracketKeyword_4()); 

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
    // $ANTLR end "rule__ElkLabel__Group__4__Impl"


    // $ANTLR start "rule__ElkLabel__Group_2__0"
    // InternalElkGraphJson.g:2861:1: rule__ElkLabel__Group_2__0 : rule__ElkLabel__Group_2__0__Impl rule__ElkLabel__Group_2__1 ;
    public final void rule__ElkLabel__Group_2__0() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:2865:1: ( rule__ElkLabel__Group_2__0__Impl rule__ElkLabel__Group_2__1 )
            // InternalElkGraphJson.g:2866:2: rule__ElkLabel__Group_2__0__Impl rule__ElkLabel__Group_2__1
            {
            pushFollow(FOLLOW_11);
            rule__ElkLabel__Group_2__0__Impl();

            state._fsp--;

            pushFollow(FOLLOW_2);
            rule__ElkLabel__Group_2__1();

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
    // $ANTLR end "rule__ElkLabel__Group_2__0"


    // $ANTLR start "rule__ElkLabel__Group_2__0__Impl"
    // InternalElkGraphJson.g:2873:1: rule__ElkLabel__Group_2__0__Impl : ( ',' ) ;
    public final void rule__ElkLabel__Group_2__0__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:2877:1: ( ( ',' ) )
            // InternalElkGraphJson.g:2878:1: ( ',' )
            {
            // InternalElkGraphJson.g:2878:1: ( ',' )
            // InternalElkGraphJson.g:2879:2: ','
            {
             before(grammarAccess.getElkLabelAccess().getCommaKeyword_2_0()); 
            match(input,59,FOLLOW_2); 
             after(grammarAccess.getElkLabelAccess().getCommaKeyword_2_0()); 

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
    // $ANTLR end "rule__ElkLabel__Group_2__0__Impl"


    // $ANTLR start "rule__ElkLabel__Group_2__1"
    // InternalElkGraphJson.g:2888:1: rule__ElkLabel__Group_2__1 : rule__ElkLabel__Group_2__1__Impl ;
    public final void rule__ElkLabel__Group_2__1() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:2892:1: ( rule__ElkLabel__Group_2__1__Impl )
            // InternalElkGraphJson.g:2893:2: rule__ElkLabel__Group_2__1__Impl
            {
            pushFollow(FOLLOW_2);
            rule__ElkLabel__Group_2__1__Impl();

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
    // $ANTLR end "rule__ElkLabel__Group_2__1"


    // $ANTLR start "rule__ElkLabel__Group_2__1__Impl"
    // InternalElkGraphJson.g:2899:1: rule__ElkLabel__Group_2__1__Impl : ( ruleLabelElement ) ;
    public final void rule__ElkLabel__Group_2__1__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:2903:1: ( ( ruleLabelElement ) )
            // InternalElkGraphJson.g:2904:1: ( ruleLabelElement )
            {
            // InternalElkGraphJson.g:2904:1: ( ruleLabelElement )
            // InternalElkGraphJson.g:2905:2: ruleLabelElement
            {
             before(grammarAccess.getElkLabelAccess().getLabelElementParserRuleCall_2_1()); 
            pushFollow(FOLLOW_2);
            ruleLabelElement();

            state._fsp--;

             after(grammarAccess.getElkLabelAccess().getLabelElementParserRuleCall_2_1()); 

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
    // $ANTLR end "rule__ElkLabel__Group_2__1__Impl"


    // $ANTLR start "rule__LabelElement__Group_2__0"
    // InternalElkGraphJson.g:2915:1: rule__LabelElement__Group_2__0 : rule__LabelElement__Group_2__0__Impl rule__LabelElement__Group_2__1 ;
    public final void rule__LabelElement__Group_2__0() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:2919:1: ( rule__LabelElement__Group_2__0__Impl rule__LabelElement__Group_2__1 )
            // InternalElkGraphJson.g:2920:2: rule__LabelElement__Group_2__0__Impl rule__LabelElement__Group_2__1
            {
            pushFollow(FOLLOW_8);
            rule__LabelElement__Group_2__0__Impl();

            state._fsp--;

            pushFollow(FOLLOW_2);
            rule__LabelElement__Group_2__1();

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
    // $ANTLR end "rule__LabelElement__Group_2__0"


    // $ANTLR start "rule__LabelElement__Group_2__0__Impl"
    // InternalElkGraphJson.g:2927:1: rule__LabelElement__Group_2__0__Impl : ( ruleKeyText ) ;
    public final void rule__LabelElement__Group_2__0__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:2931:1: ( ( ruleKeyText ) )
            // InternalElkGraphJson.g:2932:1: ( ruleKeyText )
            {
            // InternalElkGraphJson.g:2932:1: ( ruleKeyText )
            // InternalElkGraphJson.g:2933:2: ruleKeyText
            {
             before(grammarAccess.getLabelElementAccess().getKeyTextParserRuleCall_2_0()); 
            pushFollow(FOLLOW_2);
            ruleKeyText();

            state._fsp--;

             after(grammarAccess.getLabelElementAccess().getKeyTextParserRuleCall_2_0()); 

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
    // $ANTLR end "rule__LabelElement__Group_2__0__Impl"


    // $ANTLR start "rule__LabelElement__Group_2__1"
    // InternalElkGraphJson.g:2942:1: rule__LabelElement__Group_2__1 : rule__LabelElement__Group_2__1__Impl rule__LabelElement__Group_2__2 ;
    public final void rule__LabelElement__Group_2__1() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:2946:1: ( rule__LabelElement__Group_2__1__Impl rule__LabelElement__Group_2__2 )
            // InternalElkGraphJson.g:2947:2: rule__LabelElement__Group_2__1__Impl rule__LabelElement__Group_2__2
            {
            pushFollow(FOLLOW_12);
            rule__LabelElement__Group_2__1__Impl();

            state._fsp--;

            pushFollow(FOLLOW_2);
            rule__LabelElement__Group_2__2();

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
    // $ANTLR end "rule__LabelElement__Group_2__1"


    // $ANTLR start "rule__LabelElement__Group_2__1__Impl"
    // InternalElkGraphJson.g:2954:1: rule__LabelElement__Group_2__1__Impl : ( ':' ) ;
    public final void rule__LabelElement__Group_2__1__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:2958:1: ( ( ':' ) )
            // InternalElkGraphJson.g:2959:1: ( ':' )
            {
            // InternalElkGraphJson.g:2959:1: ( ':' )
            // InternalElkGraphJson.g:2960:2: ':'
            {
             before(grammarAccess.getLabelElementAccess().getColonKeyword_2_1()); 
            match(input,61,FOLLOW_2); 
             after(grammarAccess.getLabelElementAccess().getColonKeyword_2_1()); 

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
    // $ANTLR end "rule__LabelElement__Group_2__1__Impl"


    // $ANTLR start "rule__LabelElement__Group_2__2"
    // InternalElkGraphJson.g:2969:1: rule__LabelElement__Group_2__2 : rule__LabelElement__Group_2__2__Impl ;
    public final void rule__LabelElement__Group_2__2() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:2973:1: ( rule__LabelElement__Group_2__2__Impl )
            // InternalElkGraphJson.g:2974:2: rule__LabelElement__Group_2__2__Impl
            {
            pushFollow(FOLLOW_2);
            rule__LabelElement__Group_2__2__Impl();

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
    // $ANTLR end "rule__LabelElement__Group_2__2"


    // $ANTLR start "rule__LabelElement__Group_2__2__Impl"
    // InternalElkGraphJson.g:2980:1: rule__LabelElement__Group_2__2__Impl : ( ( rule__LabelElement__TextAssignment_2_2 ) ) ;
    public final void rule__LabelElement__Group_2__2__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:2984:1: ( ( ( rule__LabelElement__TextAssignment_2_2 ) ) )
            // InternalElkGraphJson.g:2985:1: ( ( rule__LabelElement__TextAssignment_2_2 ) )
            {
            // InternalElkGraphJson.g:2985:1: ( ( rule__LabelElement__TextAssignment_2_2 ) )
            // InternalElkGraphJson.g:2986:2: ( rule__LabelElement__TextAssignment_2_2 )
            {
             before(grammarAccess.getLabelElementAccess().getTextAssignment_2_2()); 
            // InternalElkGraphJson.g:2987:2: ( rule__LabelElement__TextAssignment_2_2 )
            // InternalElkGraphJson.g:2987:3: rule__LabelElement__TextAssignment_2_2
            {
            pushFollow(FOLLOW_2);
            rule__LabelElement__TextAssignment_2_2();

            state._fsp--;


            }

             after(grammarAccess.getLabelElementAccess().getTextAssignment_2_2()); 

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
    // $ANTLR end "rule__LabelElement__Group_2__2__Impl"


    // $ANTLR start "rule__LabelElement__Group_3__0"
    // InternalElkGraphJson.g:2996:1: rule__LabelElement__Group_3__0 : rule__LabelElement__Group_3__0__Impl rule__LabelElement__Group_3__1 ;
    public final void rule__LabelElement__Group_3__0() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:3000:1: ( rule__LabelElement__Group_3__0__Impl rule__LabelElement__Group_3__1 )
            // InternalElkGraphJson.g:3001:2: rule__LabelElement__Group_3__0__Impl rule__LabelElement__Group_3__1
            {
            pushFollow(FOLLOW_8);
            rule__LabelElement__Group_3__0__Impl();

            state._fsp--;

            pushFollow(FOLLOW_2);
            rule__LabelElement__Group_3__1();

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
    // $ANTLR end "rule__LabelElement__Group_3__0"


    // $ANTLR start "rule__LabelElement__Group_3__0__Impl"
    // InternalElkGraphJson.g:3008:1: rule__LabelElement__Group_3__0__Impl : ( ruleKeyLabels ) ;
    public final void rule__LabelElement__Group_3__0__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:3012:1: ( ( ruleKeyLabels ) )
            // InternalElkGraphJson.g:3013:1: ( ruleKeyLabels )
            {
            // InternalElkGraphJson.g:3013:1: ( ruleKeyLabels )
            // InternalElkGraphJson.g:3014:2: ruleKeyLabels
            {
             before(grammarAccess.getLabelElementAccess().getKeyLabelsParserRuleCall_3_0()); 
            pushFollow(FOLLOW_2);
            ruleKeyLabels();

            state._fsp--;

             after(grammarAccess.getLabelElementAccess().getKeyLabelsParserRuleCall_3_0()); 

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
    // $ANTLR end "rule__LabelElement__Group_3__0__Impl"


    // $ANTLR start "rule__LabelElement__Group_3__1"
    // InternalElkGraphJson.g:3023:1: rule__LabelElement__Group_3__1 : rule__LabelElement__Group_3__1__Impl rule__LabelElement__Group_3__2 ;
    public final void rule__LabelElement__Group_3__1() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:3027:1: ( rule__LabelElement__Group_3__1__Impl rule__LabelElement__Group_3__2 )
            // InternalElkGraphJson.g:3028:2: rule__LabelElement__Group_3__1__Impl rule__LabelElement__Group_3__2
            {
            pushFollow(FOLLOW_9);
            rule__LabelElement__Group_3__1__Impl();

            state._fsp--;

            pushFollow(FOLLOW_2);
            rule__LabelElement__Group_3__2();

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
    // $ANTLR end "rule__LabelElement__Group_3__1"


    // $ANTLR start "rule__LabelElement__Group_3__1__Impl"
    // InternalElkGraphJson.g:3035:1: rule__LabelElement__Group_3__1__Impl : ( ':' ) ;
    public final void rule__LabelElement__Group_3__1__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:3039:1: ( ( ':' ) )
            // InternalElkGraphJson.g:3040:1: ( ':' )
            {
            // InternalElkGraphJson.g:3040:1: ( ':' )
            // InternalElkGraphJson.g:3041:2: ':'
            {
             before(grammarAccess.getLabelElementAccess().getColonKeyword_3_1()); 
            match(input,61,FOLLOW_2); 
             after(grammarAccess.getLabelElementAccess().getColonKeyword_3_1()); 

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
    // $ANTLR end "rule__LabelElement__Group_3__1__Impl"


    // $ANTLR start "rule__LabelElement__Group_3__2"
    // InternalElkGraphJson.g:3050:1: rule__LabelElement__Group_3__2 : rule__LabelElement__Group_3__2__Impl ;
    public final void rule__LabelElement__Group_3__2() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:3054:1: ( rule__LabelElement__Group_3__2__Impl )
            // InternalElkGraphJson.g:3055:2: rule__LabelElement__Group_3__2__Impl
            {
            pushFollow(FOLLOW_2);
            rule__LabelElement__Group_3__2__Impl();

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
    // $ANTLR end "rule__LabelElement__Group_3__2"


    // $ANTLR start "rule__LabelElement__Group_3__2__Impl"
    // InternalElkGraphJson.g:3061:1: rule__LabelElement__Group_3__2__Impl : ( ruleElkGraphElementLabels ) ;
    public final void rule__LabelElement__Group_3__2__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:3065:1: ( ( ruleElkGraphElementLabels ) )
            // InternalElkGraphJson.g:3066:1: ( ruleElkGraphElementLabels )
            {
            // InternalElkGraphJson.g:3066:1: ( ruleElkGraphElementLabels )
            // InternalElkGraphJson.g:3067:2: ruleElkGraphElementLabels
            {
             before(grammarAccess.getLabelElementAccess().getElkGraphElementLabelsParserRuleCall_3_2()); 
            pushFollow(FOLLOW_2);
            ruleElkGraphElementLabels();

            state._fsp--;

             after(grammarAccess.getLabelElementAccess().getElkGraphElementLabelsParserRuleCall_3_2()); 

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
    // $ANTLR end "rule__LabelElement__Group_3__2__Impl"


    // $ANTLR start "rule__LabelElement__Group_4__0"
    // InternalElkGraphJson.g:3077:1: rule__LabelElement__Group_4__0 : rule__LabelElement__Group_4__0__Impl rule__LabelElement__Group_4__1 ;
    public final void rule__LabelElement__Group_4__0() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:3081:1: ( rule__LabelElement__Group_4__0__Impl rule__LabelElement__Group_4__1 )
            // InternalElkGraphJson.g:3082:2: rule__LabelElement__Group_4__0__Impl rule__LabelElement__Group_4__1
            {
            pushFollow(FOLLOW_8);
            rule__LabelElement__Group_4__0__Impl();

            state._fsp--;

            pushFollow(FOLLOW_2);
            rule__LabelElement__Group_4__1();

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
    // $ANTLR end "rule__LabelElement__Group_4__0"


    // $ANTLR start "rule__LabelElement__Group_4__0__Impl"
    // InternalElkGraphJson.g:3089:1: rule__LabelElement__Group_4__0__Impl : ( ruleKeyLayoutOptions ) ;
    public final void rule__LabelElement__Group_4__0__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:3093:1: ( ( ruleKeyLayoutOptions ) )
            // InternalElkGraphJson.g:3094:1: ( ruleKeyLayoutOptions )
            {
            // InternalElkGraphJson.g:3094:1: ( ruleKeyLayoutOptions )
            // InternalElkGraphJson.g:3095:2: ruleKeyLayoutOptions
            {
             before(grammarAccess.getLabelElementAccess().getKeyLayoutOptionsParserRuleCall_4_0()); 
            pushFollow(FOLLOW_2);
            ruleKeyLayoutOptions();

            state._fsp--;

             after(grammarAccess.getLabelElementAccess().getKeyLayoutOptionsParserRuleCall_4_0()); 

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
    // $ANTLR end "rule__LabelElement__Group_4__0__Impl"


    // $ANTLR start "rule__LabelElement__Group_4__1"
    // InternalElkGraphJson.g:3104:1: rule__LabelElement__Group_4__1 : rule__LabelElement__Group_4__1__Impl rule__LabelElement__Group_4__2 ;
    public final void rule__LabelElement__Group_4__1() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:3108:1: ( rule__LabelElement__Group_4__1__Impl rule__LabelElement__Group_4__2 )
            // InternalElkGraphJson.g:3109:2: rule__LabelElement__Group_4__1__Impl rule__LabelElement__Group_4__2
            {
            pushFollow(FOLLOW_3);
            rule__LabelElement__Group_4__1__Impl();

            state._fsp--;

            pushFollow(FOLLOW_2);
            rule__LabelElement__Group_4__2();

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
    // $ANTLR end "rule__LabelElement__Group_4__1"


    // $ANTLR start "rule__LabelElement__Group_4__1__Impl"
    // InternalElkGraphJson.g:3116:1: rule__LabelElement__Group_4__1__Impl : ( ':' ) ;
    public final void rule__LabelElement__Group_4__1__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:3120:1: ( ( ':' ) )
            // InternalElkGraphJson.g:3121:1: ( ':' )
            {
            // InternalElkGraphJson.g:3121:1: ( ':' )
            // InternalElkGraphJson.g:3122:2: ':'
            {
             before(grammarAccess.getLabelElementAccess().getColonKeyword_4_1()); 
            match(input,61,FOLLOW_2); 
             after(grammarAccess.getLabelElementAccess().getColonKeyword_4_1()); 

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
    // $ANTLR end "rule__LabelElement__Group_4__1__Impl"


    // $ANTLR start "rule__LabelElement__Group_4__2"
    // InternalElkGraphJson.g:3131:1: rule__LabelElement__Group_4__2 : rule__LabelElement__Group_4__2__Impl ;
    public final void rule__LabelElement__Group_4__2() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:3135:1: ( rule__LabelElement__Group_4__2__Impl )
            // InternalElkGraphJson.g:3136:2: rule__LabelElement__Group_4__2__Impl
            {
            pushFollow(FOLLOW_2);
            rule__LabelElement__Group_4__2__Impl();

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
    // $ANTLR end "rule__LabelElement__Group_4__2"


    // $ANTLR start "rule__LabelElement__Group_4__2__Impl"
    // InternalElkGraphJson.g:3142:1: rule__LabelElement__Group_4__2__Impl : ( ruleElkGraphElementProperties ) ;
    public final void rule__LabelElement__Group_4__2__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:3146:1: ( ( ruleElkGraphElementProperties ) )
            // InternalElkGraphJson.g:3147:1: ( ruleElkGraphElementProperties )
            {
            // InternalElkGraphJson.g:3147:1: ( ruleElkGraphElementProperties )
            // InternalElkGraphJson.g:3148:2: ruleElkGraphElementProperties
            {
             before(grammarAccess.getLabelElementAccess().getElkGraphElementPropertiesParserRuleCall_4_2()); 
            pushFollow(FOLLOW_2);
            ruleElkGraphElementProperties();

            state._fsp--;

             after(grammarAccess.getLabelElementAccess().getElkGraphElementPropertiesParserRuleCall_4_2()); 

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
    // $ANTLR end "rule__LabelElement__Group_4__2__Impl"


    // $ANTLR start "rule__ElkEdge__Group__0"
    // InternalElkGraphJson.g:3158:1: rule__ElkEdge__Group__0 : rule__ElkEdge__Group__0__Impl rule__ElkEdge__Group__1 ;
    public final void rule__ElkEdge__Group__0() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:3162:1: ( rule__ElkEdge__Group__0__Impl rule__ElkEdge__Group__1 )
            // InternalElkGraphJson.g:3163:2: rule__ElkEdge__Group__0__Impl rule__ElkEdge__Group__1
            {
            pushFollow(FOLLOW_13);
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
    // InternalElkGraphJson.g:3170:1: rule__ElkEdge__Group__0__Impl : ( '{' ) ;
    public final void rule__ElkEdge__Group__0__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:3174:1: ( ( '{' ) )
            // InternalElkGraphJson.g:3175:1: ( '{' )
            {
            // InternalElkGraphJson.g:3175:1: ( '{' )
            // InternalElkGraphJson.g:3176:2: '{'
            {
             before(grammarAccess.getElkEdgeAccess().getLeftCurlyBracketKeyword_0()); 
            match(input,58,FOLLOW_2); 
             after(grammarAccess.getElkEdgeAccess().getLeftCurlyBracketKeyword_0()); 

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
    // InternalElkGraphJson.g:3185:1: rule__ElkEdge__Group__1 : rule__ElkEdge__Group__1__Impl rule__ElkEdge__Group__2 ;
    public final void rule__ElkEdge__Group__1() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:3189:1: ( rule__ElkEdge__Group__1__Impl rule__ElkEdge__Group__2 )
            // InternalElkGraphJson.g:3190:2: rule__ElkEdge__Group__1__Impl rule__ElkEdge__Group__2
            {
            pushFollow(FOLLOW_10);
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
    // InternalElkGraphJson.g:3197:1: rule__ElkEdge__Group__1__Impl : ( ruleEdgeElement ) ;
    public final void rule__ElkEdge__Group__1__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:3201:1: ( ( ruleEdgeElement ) )
            // InternalElkGraphJson.g:3202:1: ( ruleEdgeElement )
            {
            // InternalElkGraphJson.g:3202:1: ( ruleEdgeElement )
            // InternalElkGraphJson.g:3203:2: ruleEdgeElement
            {
             before(grammarAccess.getElkEdgeAccess().getEdgeElementParserRuleCall_1()); 
            pushFollow(FOLLOW_2);
            ruleEdgeElement();

            state._fsp--;

             after(grammarAccess.getElkEdgeAccess().getEdgeElementParserRuleCall_1()); 

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
    // InternalElkGraphJson.g:3212:1: rule__ElkEdge__Group__2 : rule__ElkEdge__Group__2__Impl rule__ElkEdge__Group__3 ;
    public final void rule__ElkEdge__Group__2() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:3216:1: ( rule__ElkEdge__Group__2__Impl rule__ElkEdge__Group__3 )
            // InternalElkGraphJson.g:3217:2: rule__ElkEdge__Group__2__Impl rule__ElkEdge__Group__3
            {
            pushFollow(FOLLOW_10);
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
    // InternalElkGraphJson.g:3224:1: rule__ElkEdge__Group__2__Impl : ( ( rule__ElkEdge__Group_2__0 )* ) ;
    public final void rule__ElkEdge__Group__2__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:3228:1: ( ( ( rule__ElkEdge__Group_2__0 )* ) )
            // InternalElkGraphJson.g:3229:1: ( ( rule__ElkEdge__Group_2__0 )* )
            {
            // InternalElkGraphJson.g:3229:1: ( ( rule__ElkEdge__Group_2__0 )* )
            // InternalElkGraphJson.g:3230:2: ( rule__ElkEdge__Group_2__0 )*
            {
             before(grammarAccess.getElkEdgeAccess().getGroup_2()); 
            // InternalElkGraphJson.g:3231:2: ( rule__ElkEdge__Group_2__0 )*
            loop33:
            do {
                int alt33=2;
                int LA33_0 = input.LA(1);

                if ( (LA33_0==59) ) {
                    int LA33_1 = input.LA(2);

                    if ( ((LA33_1>=RULE_STRING && LA33_1<=RULE_ID)||(LA33_1>=22 && LA33_1<=24)||(LA33_1>=28 && LA33_1<=36)||(LA33_1>=49 && LA33_1<=54)) ) {
                        alt33=1;
                    }


                }


                switch (alt33) {
            	case 1 :
            	    // InternalElkGraphJson.g:3231:3: rule__ElkEdge__Group_2__0
            	    {
            	    pushFollow(FOLLOW_6);
            	    rule__ElkEdge__Group_2__0();

            	    state._fsp--;


            	    }
            	    break;

            	default :
            	    break loop33;
                }
            } while (true);

             after(grammarAccess.getElkEdgeAccess().getGroup_2()); 

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
    // InternalElkGraphJson.g:3239:1: rule__ElkEdge__Group__3 : rule__ElkEdge__Group__3__Impl rule__ElkEdge__Group__4 ;
    public final void rule__ElkEdge__Group__3() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:3243:1: ( rule__ElkEdge__Group__3__Impl rule__ElkEdge__Group__4 )
            // InternalElkGraphJson.g:3244:2: rule__ElkEdge__Group__3__Impl rule__ElkEdge__Group__4
            {
            pushFollow(FOLLOW_10);
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
    // InternalElkGraphJson.g:3251:1: rule__ElkEdge__Group__3__Impl : ( ( ',' )? ) ;
    public final void rule__ElkEdge__Group__3__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:3255:1: ( ( ( ',' )? ) )
            // InternalElkGraphJson.g:3256:1: ( ( ',' )? )
            {
            // InternalElkGraphJson.g:3256:1: ( ( ',' )? )
            // InternalElkGraphJson.g:3257:2: ( ',' )?
            {
             before(grammarAccess.getElkEdgeAccess().getCommaKeyword_3()); 
            // InternalElkGraphJson.g:3258:2: ( ',' )?
            int alt34=2;
            int LA34_0 = input.LA(1);

            if ( (LA34_0==59) ) {
                alt34=1;
            }
            switch (alt34) {
                case 1 :
                    // InternalElkGraphJson.g:3258:3: ','
                    {
                    match(input,59,FOLLOW_2); 

                    }
                    break;

            }

             after(grammarAccess.getElkEdgeAccess().getCommaKeyword_3()); 

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
    // InternalElkGraphJson.g:3266:1: rule__ElkEdge__Group__4 : rule__ElkEdge__Group__4__Impl ;
    public final void rule__ElkEdge__Group__4() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:3270:1: ( rule__ElkEdge__Group__4__Impl )
            // InternalElkGraphJson.g:3271:2: rule__ElkEdge__Group__4__Impl
            {
            pushFollow(FOLLOW_2);
            rule__ElkEdge__Group__4__Impl();

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
    // InternalElkGraphJson.g:3277:1: rule__ElkEdge__Group__4__Impl : ( '}' ) ;
    public final void rule__ElkEdge__Group__4__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:3281:1: ( ( '}' ) )
            // InternalElkGraphJson.g:3282:1: ( '}' )
            {
            // InternalElkGraphJson.g:3282:1: ( '}' )
            // InternalElkGraphJson.g:3283:2: '}'
            {
             before(grammarAccess.getElkEdgeAccess().getRightCurlyBracketKeyword_4()); 
            match(input,60,FOLLOW_2); 
             after(grammarAccess.getElkEdgeAccess().getRightCurlyBracketKeyword_4()); 

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


    // $ANTLR start "rule__ElkEdge__Group_2__0"
    // InternalElkGraphJson.g:3293:1: rule__ElkEdge__Group_2__0 : rule__ElkEdge__Group_2__0__Impl rule__ElkEdge__Group_2__1 ;
    public final void rule__ElkEdge__Group_2__0() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:3297:1: ( rule__ElkEdge__Group_2__0__Impl rule__ElkEdge__Group_2__1 )
            // InternalElkGraphJson.g:3298:2: rule__ElkEdge__Group_2__0__Impl rule__ElkEdge__Group_2__1
            {
            pushFollow(FOLLOW_13);
            rule__ElkEdge__Group_2__0__Impl();

            state._fsp--;

            pushFollow(FOLLOW_2);
            rule__ElkEdge__Group_2__1();

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
    // $ANTLR end "rule__ElkEdge__Group_2__0"


    // $ANTLR start "rule__ElkEdge__Group_2__0__Impl"
    // InternalElkGraphJson.g:3305:1: rule__ElkEdge__Group_2__0__Impl : ( ',' ) ;
    public final void rule__ElkEdge__Group_2__0__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:3309:1: ( ( ',' ) )
            // InternalElkGraphJson.g:3310:1: ( ',' )
            {
            // InternalElkGraphJson.g:3310:1: ( ',' )
            // InternalElkGraphJson.g:3311:2: ','
            {
             before(grammarAccess.getElkEdgeAccess().getCommaKeyword_2_0()); 
            match(input,59,FOLLOW_2); 
             after(grammarAccess.getElkEdgeAccess().getCommaKeyword_2_0()); 

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
    // $ANTLR end "rule__ElkEdge__Group_2__0__Impl"


    // $ANTLR start "rule__ElkEdge__Group_2__1"
    // InternalElkGraphJson.g:3320:1: rule__ElkEdge__Group_2__1 : rule__ElkEdge__Group_2__1__Impl ;
    public final void rule__ElkEdge__Group_2__1() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:3324:1: ( rule__ElkEdge__Group_2__1__Impl )
            // InternalElkGraphJson.g:3325:2: rule__ElkEdge__Group_2__1__Impl
            {
            pushFollow(FOLLOW_2);
            rule__ElkEdge__Group_2__1__Impl();

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
    // $ANTLR end "rule__ElkEdge__Group_2__1"


    // $ANTLR start "rule__ElkEdge__Group_2__1__Impl"
    // InternalElkGraphJson.g:3331:1: rule__ElkEdge__Group_2__1__Impl : ( ruleEdgeElement ) ;
    public final void rule__ElkEdge__Group_2__1__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:3335:1: ( ( ruleEdgeElement ) )
            // InternalElkGraphJson.g:3336:1: ( ruleEdgeElement )
            {
            // InternalElkGraphJson.g:3336:1: ( ruleEdgeElement )
            // InternalElkGraphJson.g:3337:2: ruleEdgeElement
            {
             before(grammarAccess.getElkEdgeAccess().getEdgeElementParserRuleCall_2_1()); 
            pushFollow(FOLLOW_2);
            ruleEdgeElement();

            state._fsp--;

             after(grammarAccess.getElkEdgeAccess().getEdgeElementParserRuleCall_2_1()); 

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
    // $ANTLR end "rule__ElkEdge__Group_2__1__Impl"


    // $ANTLR start "rule__EdgeElement__Group_1__0"
    // InternalElkGraphJson.g:3347:1: rule__EdgeElement__Group_1__0 : rule__EdgeElement__Group_1__0__Impl rule__EdgeElement__Group_1__1 ;
    public final void rule__EdgeElement__Group_1__0() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:3351:1: ( rule__EdgeElement__Group_1__0__Impl rule__EdgeElement__Group_1__1 )
            // InternalElkGraphJson.g:3352:2: rule__EdgeElement__Group_1__0__Impl rule__EdgeElement__Group_1__1
            {
            pushFollow(FOLLOW_8);
            rule__EdgeElement__Group_1__0__Impl();

            state._fsp--;

            pushFollow(FOLLOW_2);
            rule__EdgeElement__Group_1__1();

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
    // $ANTLR end "rule__EdgeElement__Group_1__0"


    // $ANTLR start "rule__EdgeElement__Group_1__0__Impl"
    // InternalElkGraphJson.g:3359:1: rule__EdgeElement__Group_1__0__Impl : ( ruleKeySources ) ;
    public final void rule__EdgeElement__Group_1__0__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:3363:1: ( ( ruleKeySources ) )
            // InternalElkGraphJson.g:3364:1: ( ruleKeySources )
            {
            // InternalElkGraphJson.g:3364:1: ( ruleKeySources )
            // InternalElkGraphJson.g:3365:2: ruleKeySources
            {
             before(grammarAccess.getEdgeElementAccess().getKeySourcesParserRuleCall_1_0()); 
            pushFollow(FOLLOW_2);
            ruleKeySources();

            state._fsp--;

             after(grammarAccess.getEdgeElementAccess().getKeySourcesParserRuleCall_1_0()); 

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
    // $ANTLR end "rule__EdgeElement__Group_1__0__Impl"


    // $ANTLR start "rule__EdgeElement__Group_1__1"
    // InternalElkGraphJson.g:3374:1: rule__EdgeElement__Group_1__1 : rule__EdgeElement__Group_1__1__Impl rule__EdgeElement__Group_1__2 ;
    public final void rule__EdgeElement__Group_1__1() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:3378:1: ( rule__EdgeElement__Group_1__1__Impl rule__EdgeElement__Group_1__2 )
            // InternalElkGraphJson.g:3379:2: rule__EdgeElement__Group_1__1__Impl rule__EdgeElement__Group_1__2
            {
            pushFollow(FOLLOW_9);
            rule__EdgeElement__Group_1__1__Impl();

            state._fsp--;

            pushFollow(FOLLOW_2);
            rule__EdgeElement__Group_1__2();

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
    // $ANTLR end "rule__EdgeElement__Group_1__1"


    // $ANTLR start "rule__EdgeElement__Group_1__1__Impl"
    // InternalElkGraphJson.g:3386:1: rule__EdgeElement__Group_1__1__Impl : ( ':' ) ;
    public final void rule__EdgeElement__Group_1__1__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:3390:1: ( ( ':' ) )
            // InternalElkGraphJson.g:3391:1: ( ':' )
            {
            // InternalElkGraphJson.g:3391:1: ( ':' )
            // InternalElkGraphJson.g:3392:2: ':'
            {
             before(grammarAccess.getEdgeElementAccess().getColonKeyword_1_1()); 
            match(input,61,FOLLOW_2); 
             after(grammarAccess.getEdgeElementAccess().getColonKeyword_1_1()); 

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
    // $ANTLR end "rule__EdgeElement__Group_1__1__Impl"


    // $ANTLR start "rule__EdgeElement__Group_1__2"
    // InternalElkGraphJson.g:3401:1: rule__EdgeElement__Group_1__2 : rule__EdgeElement__Group_1__2__Impl ;
    public final void rule__EdgeElement__Group_1__2() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:3405:1: ( rule__EdgeElement__Group_1__2__Impl )
            // InternalElkGraphJson.g:3406:2: rule__EdgeElement__Group_1__2__Impl
            {
            pushFollow(FOLLOW_2);
            rule__EdgeElement__Group_1__2__Impl();

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
    // $ANTLR end "rule__EdgeElement__Group_1__2"


    // $ANTLR start "rule__EdgeElement__Group_1__2__Impl"
    // InternalElkGraphJson.g:3412:1: rule__EdgeElement__Group_1__2__Impl : ( ruleElkEdgeSources ) ;
    public final void rule__EdgeElement__Group_1__2__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:3416:1: ( ( ruleElkEdgeSources ) )
            // InternalElkGraphJson.g:3417:1: ( ruleElkEdgeSources )
            {
            // InternalElkGraphJson.g:3417:1: ( ruleElkEdgeSources )
            // InternalElkGraphJson.g:3418:2: ruleElkEdgeSources
            {
             before(grammarAccess.getEdgeElementAccess().getElkEdgeSourcesParserRuleCall_1_2()); 
            pushFollow(FOLLOW_2);
            ruleElkEdgeSources();

            state._fsp--;

             after(grammarAccess.getEdgeElementAccess().getElkEdgeSourcesParserRuleCall_1_2()); 

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
    // $ANTLR end "rule__EdgeElement__Group_1__2__Impl"


    // $ANTLR start "rule__EdgeElement__Group_2__0"
    // InternalElkGraphJson.g:3428:1: rule__EdgeElement__Group_2__0 : rule__EdgeElement__Group_2__0__Impl rule__EdgeElement__Group_2__1 ;
    public final void rule__EdgeElement__Group_2__0() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:3432:1: ( rule__EdgeElement__Group_2__0__Impl rule__EdgeElement__Group_2__1 )
            // InternalElkGraphJson.g:3433:2: rule__EdgeElement__Group_2__0__Impl rule__EdgeElement__Group_2__1
            {
            pushFollow(FOLLOW_8);
            rule__EdgeElement__Group_2__0__Impl();

            state._fsp--;

            pushFollow(FOLLOW_2);
            rule__EdgeElement__Group_2__1();

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
    // $ANTLR end "rule__EdgeElement__Group_2__0"


    // $ANTLR start "rule__EdgeElement__Group_2__0__Impl"
    // InternalElkGraphJson.g:3440:1: rule__EdgeElement__Group_2__0__Impl : ( ruleKeyTargets ) ;
    public final void rule__EdgeElement__Group_2__0__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:3444:1: ( ( ruleKeyTargets ) )
            // InternalElkGraphJson.g:3445:1: ( ruleKeyTargets )
            {
            // InternalElkGraphJson.g:3445:1: ( ruleKeyTargets )
            // InternalElkGraphJson.g:3446:2: ruleKeyTargets
            {
             before(grammarAccess.getEdgeElementAccess().getKeyTargetsParserRuleCall_2_0()); 
            pushFollow(FOLLOW_2);
            ruleKeyTargets();

            state._fsp--;

             after(grammarAccess.getEdgeElementAccess().getKeyTargetsParserRuleCall_2_0()); 

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
    // $ANTLR end "rule__EdgeElement__Group_2__0__Impl"


    // $ANTLR start "rule__EdgeElement__Group_2__1"
    // InternalElkGraphJson.g:3455:1: rule__EdgeElement__Group_2__1 : rule__EdgeElement__Group_2__1__Impl rule__EdgeElement__Group_2__2 ;
    public final void rule__EdgeElement__Group_2__1() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:3459:1: ( rule__EdgeElement__Group_2__1__Impl rule__EdgeElement__Group_2__2 )
            // InternalElkGraphJson.g:3460:2: rule__EdgeElement__Group_2__1__Impl rule__EdgeElement__Group_2__2
            {
            pushFollow(FOLLOW_9);
            rule__EdgeElement__Group_2__1__Impl();

            state._fsp--;

            pushFollow(FOLLOW_2);
            rule__EdgeElement__Group_2__2();

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
    // $ANTLR end "rule__EdgeElement__Group_2__1"


    // $ANTLR start "rule__EdgeElement__Group_2__1__Impl"
    // InternalElkGraphJson.g:3467:1: rule__EdgeElement__Group_2__1__Impl : ( ':' ) ;
    public final void rule__EdgeElement__Group_2__1__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:3471:1: ( ( ':' ) )
            // InternalElkGraphJson.g:3472:1: ( ':' )
            {
            // InternalElkGraphJson.g:3472:1: ( ':' )
            // InternalElkGraphJson.g:3473:2: ':'
            {
             before(grammarAccess.getEdgeElementAccess().getColonKeyword_2_1()); 
            match(input,61,FOLLOW_2); 
             after(grammarAccess.getEdgeElementAccess().getColonKeyword_2_1()); 

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
    // $ANTLR end "rule__EdgeElement__Group_2__1__Impl"


    // $ANTLR start "rule__EdgeElement__Group_2__2"
    // InternalElkGraphJson.g:3482:1: rule__EdgeElement__Group_2__2 : rule__EdgeElement__Group_2__2__Impl ;
    public final void rule__EdgeElement__Group_2__2() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:3486:1: ( rule__EdgeElement__Group_2__2__Impl )
            // InternalElkGraphJson.g:3487:2: rule__EdgeElement__Group_2__2__Impl
            {
            pushFollow(FOLLOW_2);
            rule__EdgeElement__Group_2__2__Impl();

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
    // $ANTLR end "rule__EdgeElement__Group_2__2"


    // $ANTLR start "rule__EdgeElement__Group_2__2__Impl"
    // InternalElkGraphJson.g:3493:1: rule__EdgeElement__Group_2__2__Impl : ( ruleElkEdgeTargets ) ;
    public final void rule__EdgeElement__Group_2__2__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:3497:1: ( ( ruleElkEdgeTargets ) )
            // InternalElkGraphJson.g:3498:1: ( ruleElkEdgeTargets )
            {
            // InternalElkGraphJson.g:3498:1: ( ruleElkEdgeTargets )
            // InternalElkGraphJson.g:3499:2: ruleElkEdgeTargets
            {
             before(grammarAccess.getEdgeElementAccess().getElkEdgeTargetsParserRuleCall_2_2()); 
            pushFollow(FOLLOW_2);
            ruleElkEdgeTargets();

            state._fsp--;

             after(grammarAccess.getEdgeElementAccess().getElkEdgeTargetsParserRuleCall_2_2()); 

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
    // $ANTLR end "rule__EdgeElement__Group_2__2__Impl"


    // $ANTLR start "rule__EdgeElement__Group_3__0"
    // InternalElkGraphJson.g:3509:1: rule__EdgeElement__Group_3__0 : rule__EdgeElement__Group_3__0__Impl rule__EdgeElement__Group_3__1 ;
    public final void rule__EdgeElement__Group_3__0() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:3513:1: ( rule__EdgeElement__Group_3__0__Impl rule__EdgeElement__Group_3__1 )
            // InternalElkGraphJson.g:3514:2: rule__EdgeElement__Group_3__0__Impl rule__EdgeElement__Group_3__1
            {
            pushFollow(FOLLOW_8);
            rule__EdgeElement__Group_3__0__Impl();

            state._fsp--;

            pushFollow(FOLLOW_2);
            rule__EdgeElement__Group_3__1();

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
    // $ANTLR end "rule__EdgeElement__Group_3__0"


    // $ANTLR start "rule__EdgeElement__Group_3__0__Impl"
    // InternalElkGraphJson.g:3521:1: rule__EdgeElement__Group_3__0__Impl : ( ruleKeyLabels ) ;
    public final void rule__EdgeElement__Group_3__0__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:3525:1: ( ( ruleKeyLabels ) )
            // InternalElkGraphJson.g:3526:1: ( ruleKeyLabels )
            {
            // InternalElkGraphJson.g:3526:1: ( ruleKeyLabels )
            // InternalElkGraphJson.g:3527:2: ruleKeyLabels
            {
             before(grammarAccess.getEdgeElementAccess().getKeyLabelsParserRuleCall_3_0()); 
            pushFollow(FOLLOW_2);
            ruleKeyLabels();

            state._fsp--;

             after(grammarAccess.getEdgeElementAccess().getKeyLabelsParserRuleCall_3_0()); 

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
    // $ANTLR end "rule__EdgeElement__Group_3__0__Impl"


    // $ANTLR start "rule__EdgeElement__Group_3__1"
    // InternalElkGraphJson.g:3536:1: rule__EdgeElement__Group_3__1 : rule__EdgeElement__Group_3__1__Impl rule__EdgeElement__Group_3__2 ;
    public final void rule__EdgeElement__Group_3__1() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:3540:1: ( rule__EdgeElement__Group_3__1__Impl rule__EdgeElement__Group_3__2 )
            // InternalElkGraphJson.g:3541:2: rule__EdgeElement__Group_3__1__Impl rule__EdgeElement__Group_3__2
            {
            pushFollow(FOLLOW_9);
            rule__EdgeElement__Group_3__1__Impl();

            state._fsp--;

            pushFollow(FOLLOW_2);
            rule__EdgeElement__Group_3__2();

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
    // $ANTLR end "rule__EdgeElement__Group_3__1"


    // $ANTLR start "rule__EdgeElement__Group_3__1__Impl"
    // InternalElkGraphJson.g:3548:1: rule__EdgeElement__Group_3__1__Impl : ( ':' ) ;
    public final void rule__EdgeElement__Group_3__1__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:3552:1: ( ( ':' ) )
            // InternalElkGraphJson.g:3553:1: ( ':' )
            {
            // InternalElkGraphJson.g:3553:1: ( ':' )
            // InternalElkGraphJson.g:3554:2: ':'
            {
             before(grammarAccess.getEdgeElementAccess().getColonKeyword_3_1()); 
            match(input,61,FOLLOW_2); 
             after(grammarAccess.getEdgeElementAccess().getColonKeyword_3_1()); 

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
    // $ANTLR end "rule__EdgeElement__Group_3__1__Impl"


    // $ANTLR start "rule__EdgeElement__Group_3__2"
    // InternalElkGraphJson.g:3563:1: rule__EdgeElement__Group_3__2 : rule__EdgeElement__Group_3__2__Impl ;
    public final void rule__EdgeElement__Group_3__2() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:3567:1: ( rule__EdgeElement__Group_3__2__Impl )
            // InternalElkGraphJson.g:3568:2: rule__EdgeElement__Group_3__2__Impl
            {
            pushFollow(FOLLOW_2);
            rule__EdgeElement__Group_3__2__Impl();

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
    // $ANTLR end "rule__EdgeElement__Group_3__2"


    // $ANTLR start "rule__EdgeElement__Group_3__2__Impl"
    // InternalElkGraphJson.g:3574:1: rule__EdgeElement__Group_3__2__Impl : ( ruleElkGraphElementLabels ) ;
    public final void rule__EdgeElement__Group_3__2__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:3578:1: ( ( ruleElkGraphElementLabels ) )
            // InternalElkGraphJson.g:3579:1: ( ruleElkGraphElementLabels )
            {
            // InternalElkGraphJson.g:3579:1: ( ruleElkGraphElementLabels )
            // InternalElkGraphJson.g:3580:2: ruleElkGraphElementLabels
            {
             before(grammarAccess.getEdgeElementAccess().getElkGraphElementLabelsParserRuleCall_3_2()); 
            pushFollow(FOLLOW_2);
            ruleElkGraphElementLabels();

            state._fsp--;

             after(grammarAccess.getEdgeElementAccess().getElkGraphElementLabelsParserRuleCall_3_2()); 

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
    // $ANTLR end "rule__EdgeElement__Group_3__2__Impl"


    // $ANTLR start "rule__EdgeElement__Group_4__0"
    // InternalElkGraphJson.g:3590:1: rule__EdgeElement__Group_4__0 : rule__EdgeElement__Group_4__0__Impl rule__EdgeElement__Group_4__1 ;
    public final void rule__EdgeElement__Group_4__0() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:3594:1: ( rule__EdgeElement__Group_4__0__Impl rule__EdgeElement__Group_4__1 )
            // InternalElkGraphJson.g:3595:2: rule__EdgeElement__Group_4__0__Impl rule__EdgeElement__Group_4__1
            {
            pushFollow(FOLLOW_8);
            rule__EdgeElement__Group_4__0__Impl();

            state._fsp--;

            pushFollow(FOLLOW_2);
            rule__EdgeElement__Group_4__1();

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
    // $ANTLR end "rule__EdgeElement__Group_4__0"


    // $ANTLR start "rule__EdgeElement__Group_4__0__Impl"
    // InternalElkGraphJson.g:3602:1: rule__EdgeElement__Group_4__0__Impl : ( ruleKeyLayoutOptions ) ;
    public final void rule__EdgeElement__Group_4__0__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:3606:1: ( ( ruleKeyLayoutOptions ) )
            // InternalElkGraphJson.g:3607:1: ( ruleKeyLayoutOptions )
            {
            // InternalElkGraphJson.g:3607:1: ( ruleKeyLayoutOptions )
            // InternalElkGraphJson.g:3608:2: ruleKeyLayoutOptions
            {
             before(grammarAccess.getEdgeElementAccess().getKeyLayoutOptionsParserRuleCall_4_0()); 
            pushFollow(FOLLOW_2);
            ruleKeyLayoutOptions();

            state._fsp--;

             after(grammarAccess.getEdgeElementAccess().getKeyLayoutOptionsParserRuleCall_4_0()); 

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
    // $ANTLR end "rule__EdgeElement__Group_4__0__Impl"


    // $ANTLR start "rule__EdgeElement__Group_4__1"
    // InternalElkGraphJson.g:3617:1: rule__EdgeElement__Group_4__1 : rule__EdgeElement__Group_4__1__Impl rule__EdgeElement__Group_4__2 ;
    public final void rule__EdgeElement__Group_4__1() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:3621:1: ( rule__EdgeElement__Group_4__1__Impl rule__EdgeElement__Group_4__2 )
            // InternalElkGraphJson.g:3622:2: rule__EdgeElement__Group_4__1__Impl rule__EdgeElement__Group_4__2
            {
            pushFollow(FOLLOW_3);
            rule__EdgeElement__Group_4__1__Impl();

            state._fsp--;

            pushFollow(FOLLOW_2);
            rule__EdgeElement__Group_4__2();

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
    // $ANTLR end "rule__EdgeElement__Group_4__1"


    // $ANTLR start "rule__EdgeElement__Group_4__1__Impl"
    // InternalElkGraphJson.g:3629:1: rule__EdgeElement__Group_4__1__Impl : ( ':' ) ;
    public final void rule__EdgeElement__Group_4__1__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:3633:1: ( ( ':' ) )
            // InternalElkGraphJson.g:3634:1: ( ':' )
            {
            // InternalElkGraphJson.g:3634:1: ( ':' )
            // InternalElkGraphJson.g:3635:2: ':'
            {
             before(grammarAccess.getEdgeElementAccess().getColonKeyword_4_1()); 
            match(input,61,FOLLOW_2); 
             after(grammarAccess.getEdgeElementAccess().getColonKeyword_4_1()); 

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
    // $ANTLR end "rule__EdgeElement__Group_4__1__Impl"


    // $ANTLR start "rule__EdgeElement__Group_4__2"
    // InternalElkGraphJson.g:3644:1: rule__EdgeElement__Group_4__2 : rule__EdgeElement__Group_4__2__Impl ;
    public final void rule__EdgeElement__Group_4__2() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:3648:1: ( rule__EdgeElement__Group_4__2__Impl )
            // InternalElkGraphJson.g:3649:2: rule__EdgeElement__Group_4__2__Impl
            {
            pushFollow(FOLLOW_2);
            rule__EdgeElement__Group_4__2__Impl();

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
    // $ANTLR end "rule__EdgeElement__Group_4__2"


    // $ANTLR start "rule__EdgeElement__Group_4__2__Impl"
    // InternalElkGraphJson.g:3655:1: rule__EdgeElement__Group_4__2__Impl : ( ruleElkGraphElementProperties ) ;
    public final void rule__EdgeElement__Group_4__2__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:3659:1: ( ( ruleElkGraphElementProperties ) )
            // InternalElkGraphJson.g:3660:1: ( ruleElkGraphElementProperties )
            {
            // InternalElkGraphJson.g:3660:1: ( ruleElkGraphElementProperties )
            // InternalElkGraphJson.g:3661:2: ruleElkGraphElementProperties
            {
             before(grammarAccess.getEdgeElementAccess().getElkGraphElementPropertiesParserRuleCall_4_2()); 
            pushFollow(FOLLOW_2);
            ruleElkGraphElementProperties();

            state._fsp--;

             after(grammarAccess.getEdgeElementAccess().getElkGraphElementPropertiesParserRuleCall_4_2()); 

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
    // $ANTLR end "rule__EdgeElement__Group_4__2__Impl"


    // $ANTLR start "rule__ElkEdgeSources__Group__0"
    // InternalElkGraphJson.g:3671:1: rule__ElkEdgeSources__Group__0 : rule__ElkEdgeSources__Group__0__Impl rule__ElkEdgeSources__Group__1 ;
    public final void rule__ElkEdgeSources__Group__0() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:3675:1: ( rule__ElkEdgeSources__Group__0__Impl rule__ElkEdgeSources__Group__1 )
            // InternalElkGraphJson.g:3676:2: rule__ElkEdgeSources__Group__0__Impl rule__ElkEdgeSources__Group__1
            {
            pushFollow(FOLLOW_14);
            rule__ElkEdgeSources__Group__0__Impl();

            state._fsp--;

            pushFollow(FOLLOW_2);
            rule__ElkEdgeSources__Group__1();

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
    // $ANTLR end "rule__ElkEdgeSources__Group__0"


    // $ANTLR start "rule__ElkEdgeSources__Group__0__Impl"
    // InternalElkGraphJson.g:3683:1: rule__ElkEdgeSources__Group__0__Impl : ( '[' ) ;
    public final void rule__ElkEdgeSources__Group__0__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:3687:1: ( ( '[' ) )
            // InternalElkGraphJson.g:3688:1: ( '[' )
            {
            // InternalElkGraphJson.g:3688:1: ( '[' )
            // InternalElkGraphJson.g:3689:2: '['
            {
             before(grammarAccess.getElkEdgeSourcesAccess().getLeftSquareBracketKeyword_0()); 
            match(input,62,FOLLOW_2); 
             after(grammarAccess.getElkEdgeSourcesAccess().getLeftSquareBracketKeyword_0()); 

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
    // $ANTLR end "rule__ElkEdgeSources__Group__0__Impl"


    // $ANTLR start "rule__ElkEdgeSources__Group__1"
    // InternalElkGraphJson.g:3698:1: rule__ElkEdgeSources__Group__1 : rule__ElkEdgeSources__Group__1__Impl rule__ElkEdgeSources__Group__2 ;
    public final void rule__ElkEdgeSources__Group__1() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:3702:1: ( rule__ElkEdgeSources__Group__1__Impl rule__ElkEdgeSources__Group__2 )
            // InternalElkGraphJson.g:3703:2: rule__ElkEdgeSources__Group__1__Impl rule__ElkEdgeSources__Group__2
            {
            pushFollow(FOLLOW_14);
            rule__ElkEdgeSources__Group__1__Impl();

            state._fsp--;

            pushFollow(FOLLOW_2);
            rule__ElkEdgeSources__Group__2();

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
    // $ANTLR end "rule__ElkEdgeSources__Group__1"


    // $ANTLR start "rule__ElkEdgeSources__Group__1__Impl"
    // InternalElkGraphJson.g:3710:1: rule__ElkEdgeSources__Group__1__Impl : ( ( rule__ElkEdgeSources__Group_1__0 )? ) ;
    public final void rule__ElkEdgeSources__Group__1__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:3714:1: ( ( ( rule__ElkEdgeSources__Group_1__0 )? ) )
            // InternalElkGraphJson.g:3715:1: ( ( rule__ElkEdgeSources__Group_1__0 )? )
            {
            // InternalElkGraphJson.g:3715:1: ( ( rule__ElkEdgeSources__Group_1__0 )? )
            // InternalElkGraphJson.g:3716:2: ( rule__ElkEdgeSources__Group_1__0 )?
            {
             before(grammarAccess.getElkEdgeSourcesAccess().getGroup_1()); 
            // InternalElkGraphJson.g:3717:2: ( rule__ElkEdgeSources__Group_1__0 )?
            int alt35=2;
            int LA35_0 = input.LA(1);

            if ( (LA35_0==RULE_STRING) ) {
                alt35=1;
            }
            switch (alt35) {
                case 1 :
                    // InternalElkGraphJson.g:3717:3: rule__ElkEdgeSources__Group_1__0
                    {
                    pushFollow(FOLLOW_2);
                    rule__ElkEdgeSources__Group_1__0();

                    state._fsp--;


                    }
                    break;

            }

             after(grammarAccess.getElkEdgeSourcesAccess().getGroup_1()); 

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
    // $ANTLR end "rule__ElkEdgeSources__Group__1__Impl"


    // $ANTLR start "rule__ElkEdgeSources__Group__2"
    // InternalElkGraphJson.g:3725:1: rule__ElkEdgeSources__Group__2 : rule__ElkEdgeSources__Group__2__Impl rule__ElkEdgeSources__Group__3 ;
    public final void rule__ElkEdgeSources__Group__2() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:3729:1: ( rule__ElkEdgeSources__Group__2__Impl rule__ElkEdgeSources__Group__3 )
            // InternalElkGraphJson.g:3730:2: rule__ElkEdgeSources__Group__2__Impl rule__ElkEdgeSources__Group__3
            {
            pushFollow(FOLLOW_14);
            rule__ElkEdgeSources__Group__2__Impl();

            state._fsp--;

            pushFollow(FOLLOW_2);
            rule__ElkEdgeSources__Group__3();

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
    // $ANTLR end "rule__ElkEdgeSources__Group__2"


    // $ANTLR start "rule__ElkEdgeSources__Group__2__Impl"
    // InternalElkGraphJson.g:3737:1: rule__ElkEdgeSources__Group__2__Impl : ( ( ',' )? ) ;
    public final void rule__ElkEdgeSources__Group__2__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:3741:1: ( ( ( ',' )? ) )
            // InternalElkGraphJson.g:3742:1: ( ( ',' )? )
            {
            // InternalElkGraphJson.g:3742:1: ( ( ',' )? )
            // InternalElkGraphJson.g:3743:2: ( ',' )?
            {
             before(grammarAccess.getElkEdgeSourcesAccess().getCommaKeyword_2()); 
            // InternalElkGraphJson.g:3744:2: ( ',' )?
            int alt36=2;
            int LA36_0 = input.LA(1);

            if ( (LA36_0==59) ) {
                alt36=1;
            }
            switch (alt36) {
                case 1 :
                    // InternalElkGraphJson.g:3744:3: ','
                    {
                    match(input,59,FOLLOW_2); 

                    }
                    break;

            }

             after(grammarAccess.getElkEdgeSourcesAccess().getCommaKeyword_2()); 

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
    // $ANTLR end "rule__ElkEdgeSources__Group__2__Impl"


    // $ANTLR start "rule__ElkEdgeSources__Group__3"
    // InternalElkGraphJson.g:3752:1: rule__ElkEdgeSources__Group__3 : rule__ElkEdgeSources__Group__3__Impl ;
    public final void rule__ElkEdgeSources__Group__3() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:3756:1: ( rule__ElkEdgeSources__Group__3__Impl )
            // InternalElkGraphJson.g:3757:2: rule__ElkEdgeSources__Group__3__Impl
            {
            pushFollow(FOLLOW_2);
            rule__ElkEdgeSources__Group__3__Impl();

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
    // $ANTLR end "rule__ElkEdgeSources__Group__3"


    // $ANTLR start "rule__ElkEdgeSources__Group__3__Impl"
    // InternalElkGraphJson.g:3763:1: rule__ElkEdgeSources__Group__3__Impl : ( ']' ) ;
    public final void rule__ElkEdgeSources__Group__3__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:3767:1: ( ( ']' ) )
            // InternalElkGraphJson.g:3768:1: ( ']' )
            {
            // InternalElkGraphJson.g:3768:1: ( ']' )
            // InternalElkGraphJson.g:3769:2: ']'
            {
             before(grammarAccess.getElkEdgeSourcesAccess().getRightSquareBracketKeyword_3()); 
            match(input,63,FOLLOW_2); 
             after(grammarAccess.getElkEdgeSourcesAccess().getRightSquareBracketKeyword_3()); 

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
    // $ANTLR end "rule__ElkEdgeSources__Group__3__Impl"


    // $ANTLR start "rule__ElkEdgeSources__Group_1__0"
    // InternalElkGraphJson.g:3779:1: rule__ElkEdgeSources__Group_1__0 : rule__ElkEdgeSources__Group_1__0__Impl rule__ElkEdgeSources__Group_1__1 ;
    public final void rule__ElkEdgeSources__Group_1__0() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:3783:1: ( rule__ElkEdgeSources__Group_1__0__Impl rule__ElkEdgeSources__Group_1__1 )
            // InternalElkGraphJson.g:3784:2: rule__ElkEdgeSources__Group_1__0__Impl rule__ElkEdgeSources__Group_1__1
            {
            pushFollow(FOLLOW_5);
            rule__ElkEdgeSources__Group_1__0__Impl();

            state._fsp--;

            pushFollow(FOLLOW_2);
            rule__ElkEdgeSources__Group_1__1();

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
    // $ANTLR end "rule__ElkEdgeSources__Group_1__0"


    // $ANTLR start "rule__ElkEdgeSources__Group_1__0__Impl"
    // InternalElkGraphJson.g:3791:1: rule__ElkEdgeSources__Group_1__0__Impl : ( ( rule__ElkEdgeSources__SourcesAssignment_1_0 ) ) ;
    public final void rule__ElkEdgeSources__Group_1__0__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:3795:1: ( ( ( rule__ElkEdgeSources__SourcesAssignment_1_0 ) ) )
            // InternalElkGraphJson.g:3796:1: ( ( rule__ElkEdgeSources__SourcesAssignment_1_0 ) )
            {
            // InternalElkGraphJson.g:3796:1: ( ( rule__ElkEdgeSources__SourcesAssignment_1_0 ) )
            // InternalElkGraphJson.g:3797:2: ( rule__ElkEdgeSources__SourcesAssignment_1_0 )
            {
             before(grammarAccess.getElkEdgeSourcesAccess().getSourcesAssignment_1_0()); 
            // InternalElkGraphJson.g:3798:2: ( rule__ElkEdgeSources__SourcesAssignment_1_0 )
            // InternalElkGraphJson.g:3798:3: rule__ElkEdgeSources__SourcesAssignment_1_0
            {
            pushFollow(FOLLOW_2);
            rule__ElkEdgeSources__SourcesAssignment_1_0();

            state._fsp--;


            }

             after(grammarAccess.getElkEdgeSourcesAccess().getSourcesAssignment_1_0()); 

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
    // $ANTLR end "rule__ElkEdgeSources__Group_1__0__Impl"


    // $ANTLR start "rule__ElkEdgeSources__Group_1__1"
    // InternalElkGraphJson.g:3806:1: rule__ElkEdgeSources__Group_1__1 : rule__ElkEdgeSources__Group_1__1__Impl ;
    public final void rule__ElkEdgeSources__Group_1__1() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:3810:1: ( rule__ElkEdgeSources__Group_1__1__Impl )
            // InternalElkGraphJson.g:3811:2: rule__ElkEdgeSources__Group_1__1__Impl
            {
            pushFollow(FOLLOW_2);
            rule__ElkEdgeSources__Group_1__1__Impl();

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
    // $ANTLR end "rule__ElkEdgeSources__Group_1__1"


    // $ANTLR start "rule__ElkEdgeSources__Group_1__1__Impl"
    // InternalElkGraphJson.g:3817:1: rule__ElkEdgeSources__Group_1__1__Impl : ( ( rule__ElkEdgeSources__Group_1_1__0 )* ) ;
    public final void rule__ElkEdgeSources__Group_1__1__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:3821:1: ( ( ( rule__ElkEdgeSources__Group_1_1__0 )* ) )
            // InternalElkGraphJson.g:3822:1: ( ( rule__ElkEdgeSources__Group_1_1__0 )* )
            {
            // InternalElkGraphJson.g:3822:1: ( ( rule__ElkEdgeSources__Group_1_1__0 )* )
            // InternalElkGraphJson.g:3823:2: ( rule__ElkEdgeSources__Group_1_1__0 )*
            {
             before(grammarAccess.getElkEdgeSourcesAccess().getGroup_1_1()); 
            // InternalElkGraphJson.g:3824:2: ( rule__ElkEdgeSources__Group_1_1__0 )*
            loop37:
            do {
                int alt37=2;
                int LA37_0 = input.LA(1);

                if ( (LA37_0==59) ) {
                    int LA37_1 = input.LA(2);

                    if ( (LA37_1==RULE_STRING) ) {
                        alt37=1;
                    }


                }


                switch (alt37) {
            	case 1 :
            	    // InternalElkGraphJson.g:3824:3: rule__ElkEdgeSources__Group_1_1__0
            	    {
            	    pushFollow(FOLLOW_6);
            	    rule__ElkEdgeSources__Group_1_1__0();

            	    state._fsp--;


            	    }
            	    break;

            	default :
            	    break loop37;
                }
            } while (true);

             after(grammarAccess.getElkEdgeSourcesAccess().getGroup_1_1()); 

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
    // $ANTLR end "rule__ElkEdgeSources__Group_1__1__Impl"


    // $ANTLR start "rule__ElkEdgeSources__Group_1_1__0"
    // InternalElkGraphJson.g:3833:1: rule__ElkEdgeSources__Group_1_1__0 : rule__ElkEdgeSources__Group_1_1__0__Impl rule__ElkEdgeSources__Group_1_1__1 ;
    public final void rule__ElkEdgeSources__Group_1_1__0() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:3837:1: ( rule__ElkEdgeSources__Group_1_1__0__Impl rule__ElkEdgeSources__Group_1_1__1 )
            // InternalElkGraphJson.g:3838:2: rule__ElkEdgeSources__Group_1_1__0__Impl rule__ElkEdgeSources__Group_1_1__1
            {
            pushFollow(FOLLOW_12);
            rule__ElkEdgeSources__Group_1_1__0__Impl();

            state._fsp--;

            pushFollow(FOLLOW_2);
            rule__ElkEdgeSources__Group_1_1__1();

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
    // $ANTLR end "rule__ElkEdgeSources__Group_1_1__0"


    // $ANTLR start "rule__ElkEdgeSources__Group_1_1__0__Impl"
    // InternalElkGraphJson.g:3845:1: rule__ElkEdgeSources__Group_1_1__0__Impl : ( ',' ) ;
    public final void rule__ElkEdgeSources__Group_1_1__0__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:3849:1: ( ( ',' ) )
            // InternalElkGraphJson.g:3850:1: ( ',' )
            {
            // InternalElkGraphJson.g:3850:1: ( ',' )
            // InternalElkGraphJson.g:3851:2: ','
            {
             before(grammarAccess.getElkEdgeSourcesAccess().getCommaKeyword_1_1_0()); 
            match(input,59,FOLLOW_2); 
             after(grammarAccess.getElkEdgeSourcesAccess().getCommaKeyword_1_1_0()); 

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
    // $ANTLR end "rule__ElkEdgeSources__Group_1_1__0__Impl"


    // $ANTLR start "rule__ElkEdgeSources__Group_1_1__1"
    // InternalElkGraphJson.g:3860:1: rule__ElkEdgeSources__Group_1_1__1 : rule__ElkEdgeSources__Group_1_1__1__Impl ;
    public final void rule__ElkEdgeSources__Group_1_1__1() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:3864:1: ( rule__ElkEdgeSources__Group_1_1__1__Impl )
            // InternalElkGraphJson.g:3865:2: rule__ElkEdgeSources__Group_1_1__1__Impl
            {
            pushFollow(FOLLOW_2);
            rule__ElkEdgeSources__Group_1_1__1__Impl();

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
    // $ANTLR end "rule__ElkEdgeSources__Group_1_1__1"


    // $ANTLR start "rule__ElkEdgeSources__Group_1_1__1__Impl"
    // InternalElkGraphJson.g:3871:1: rule__ElkEdgeSources__Group_1_1__1__Impl : ( ( rule__ElkEdgeSources__SourcesAssignment_1_1_1 ) ) ;
    public final void rule__ElkEdgeSources__Group_1_1__1__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:3875:1: ( ( ( rule__ElkEdgeSources__SourcesAssignment_1_1_1 ) ) )
            // InternalElkGraphJson.g:3876:1: ( ( rule__ElkEdgeSources__SourcesAssignment_1_1_1 ) )
            {
            // InternalElkGraphJson.g:3876:1: ( ( rule__ElkEdgeSources__SourcesAssignment_1_1_1 ) )
            // InternalElkGraphJson.g:3877:2: ( rule__ElkEdgeSources__SourcesAssignment_1_1_1 )
            {
             before(grammarAccess.getElkEdgeSourcesAccess().getSourcesAssignment_1_1_1()); 
            // InternalElkGraphJson.g:3878:2: ( rule__ElkEdgeSources__SourcesAssignment_1_1_1 )
            // InternalElkGraphJson.g:3878:3: rule__ElkEdgeSources__SourcesAssignment_1_1_1
            {
            pushFollow(FOLLOW_2);
            rule__ElkEdgeSources__SourcesAssignment_1_1_1();

            state._fsp--;


            }

             after(grammarAccess.getElkEdgeSourcesAccess().getSourcesAssignment_1_1_1()); 

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
    // $ANTLR end "rule__ElkEdgeSources__Group_1_1__1__Impl"


    // $ANTLR start "rule__ElkEdgeTargets__Group__0"
    // InternalElkGraphJson.g:3887:1: rule__ElkEdgeTargets__Group__0 : rule__ElkEdgeTargets__Group__0__Impl rule__ElkEdgeTargets__Group__1 ;
    public final void rule__ElkEdgeTargets__Group__0() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:3891:1: ( rule__ElkEdgeTargets__Group__0__Impl rule__ElkEdgeTargets__Group__1 )
            // InternalElkGraphJson.g:3892:2: rule__ElkEdgeTargets__Group__0__Impl rule__ElkEdgeTargets__Group__1
            {
            pushFollow(FOLLOW_14);
            rule__ElkEdgeTargets__Group__0__Impl();

            state._fsp--;

            pushFollow(FOLLOW_2);
            rule__ElkEdgeTargets__Group__1();

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
    // $ANTLR end "rule__ElkEdgeTargets__Group__0"


    // $ANTLR start "rule__ElkEdgeTargets__Group__0__Impl"
    // InternalElkGraphJson.g:3899:1: rule__ElkEdgeTargets__Group__0__Impl : ( '[' ) ;
    public final void rule__ElkEdgeTargets__Group__0__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:3903:1: ( ( '[' ) )
            // InternalElkGraphJson.g:3904:1: ( '[' )
            {
            // InternalElkGraphJson.g:3904:1: ( '[' )
            // InternalElkGraphJson.g:3905:2: '['
            {
             before(grammarAccess.getElkEdgeTargetsAccess().getLeftSquareBracketKeyword_0()); 
            match(input,62,FOLLOW_2); 
             after(grammarAccess.getElkEdgeTargetsAccess().getLeftSquareBracketKeyword_0()); 

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
    // $ANTLR end "rule__ElkEdgeTargets__Group__0__Impl"


    // $ANTLR start "rule__ElkEdgeTargets__Group__1"
    // InternalElkGraphJson.g:3914:1: rule__ElkEdgeTargets__Group__1 : rule__ElkEdgeTargets__Group__1__Impl rule__ElkEdgeTargets__Group__2 ;
    public final void rule__ElkEdgeTargets__Group__1() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:3918:1: ( rule__ElkEdgeTargets__Group__1__Impl rule__ElkEdgeTargets__Group__2 )
            // InternalElkGraphJson.g:3919:2: rule__ElkEdgeTargets__Group__1__Impl rule__ElkEdgeTargets__Group__2
            {
            pushFollow(FOLLOW_14);
            rule__ElkEdgeTargets__Group__1__Impl();

            state._fsp--;

            pushFollow(FOLLOW_2);
            rule__ElkEdgeTargets__Group__2();

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
    // $ANTLR end "rule__ElkEdgeTargets__Group__1"


    // $ANTLR start "rule__ElkEdgeTargets__Group__1__Impl"
    // InternalElkGraphJson.g:3926:1: rule__ElkEdgeTargets__Group__1__Impl : ( ( rule__ElkEdgeTargets__Group_1__0 )? ) ;
    public final void rule__ElkEdgeTargets__Group__1__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:3930:1: ( ( ( rule__ElkEdgeTargets__Group_1__0 )? ) )
            // InternalElkGraphJson.g:3931:1: ( ( rule__ElkEdgeTargets__Group_1__0 )? )
            {
            // InternalElkGraphJson.g:3931:1: ( ( rule__ElkEdgeTargets__Group_1__0 )? )
            // InternalElkGraphJson.g:3932:2: ( rule__ElkEdgeTargets__Group_1__0 )?
            {
             before(grammarAccess.getElkEdgeTargetsAccess().getGroup_1()); 
            // InternalElkGraphJson.g:3933:2: ( rule__ElkEdgeTargets__Group_1__0 )?
            int alt38=2;
            int LA38_0 = input.LA(1);

            if ( (LA38_0==RULE_STRING) ) {
                alt38=1;
            }
            switch (alt38) {
                case 1 :
                    // InternalElkGraphJson.g:3933:3: rule__ElkEdgeTargets__Group_1__0
                    {
                    pushFollow(FOLLOW_2);
                    rule__ElkEdgeTargets__Group_1__0();

                    state._fsp--;


                    }
                    break;

            }

             after(grammarAccess.getElkEdgeTargetsAccess().getGroup_1()); 

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
    // $ANTLR end "rule__ElkEdgeTargets__Group__1__Impl"


    // $ANTLR start "rule__ElkEdgeTargets__Group__2"
    // InternalElkGraphJson.g:3941:1: rule__ElkEdgeTargets__Group__2 : rule__ElkEdgeTargets__Group__2__Impl rule__ElkEdgeTargets__Group__3 ;
    public final void rule__ElkEdgeTargets__Group__2() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:3945:1: ( rule__ElkEdgeTargets__Group__2__Impl rule__ElkEdgeTargets__Group__3 )
            // InternalElkGraphJson.g:3946:2: rule__ElkEdgeTargets__Group__2__Impl rule__ElkEdgeTargets__Group__3
            {
            pushFollow(FOLLOW_14);
            rule__ElkEdgeTargets__Group__2__Impl();

            state._fsp--;

            pushFollow(FOLLOW_2);
            rule__ElkEdgeTargets__Group__3();

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
    // $ANTLR end "rule__ElkEdgeTargets__Group__2"


    // $ANTLR start "rule__ElkEdgeTargets__Group__2__Impl"
    // InternalElkGraphJson.g:3953:1: rule__ElkEdgeTargets__Group__2__Impl : ( ( ',' )? ) ;
    public final void rule__ElkEdgeTargets__Group__2__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:3957:1: ( ( ( ',' )? ) )
            // InternalElkGraphJson.g:3958:1: ( ( ',' )? )
            {
            // InternalElkGraphJson.g:3958:1: ( ( ',' )? )
            // InternalElkGraphJson.g:3959:2: ( ',' )?
            {
             before(grammarAccess.getElkEdgeTargetsAccess().getCommaKeyword_2()); 
            // InternalElkGraphJson.g:3960:2: ( ',' )?
            int alt39=2;
            int LA39_0 = input.LA(1);

            if ( (LA39_0==59) ) {
                alt39=1;
            }
            switch (alt39) {
                case 1 :
                    // InternalElkGraphJson.g:3960:3: ','
                    {
                    match(input,59,FOLLOW_2); 

                    }
                    break;

            }

             after(grammarAccess.getElkEdgeTargetsAccess().getCommaKeyword_2()); 

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
    // $ANTLR end "rule__ElkEdgeTargets__Group__2__Impl"


    // $ANTLR start "rule__ElkEdgeTargets__Group__3"
    // InternalElkGraphJson.g:3968:1: rule__ElkEdgeTargets__Group__3 : rule__ElkEdgeTargets__Group__3__Impl ;
    public final void rule__ElkEdgeTargets__Group__3() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:3972:1: ( rule__ElkEdgeTargets__Group__3__Impl )
            // InternalElkGraphJson.g:3973:2: rule__ElkEdgeTargets__Group__3__Impl
            {
            pushFollow(FOLLOW_2);
            rule__ElkEdgeTargets__Group__3__Impl();

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
    // $ANTLR end "rule__ElkEdgeTargets__Group__3"


    // $ANTLR start "rule__ElkEdgeTargets__Group__3__Impl"
    // InternalElkGraphJson.g:3979:1: rule__ElkEdgeTargets__Group__3__Impl : ( ']' ) ;
    public final void rule__ElkEdgeTargets__Group__3__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:3983:1: ( ( ']' ) )
            // InternalElkGraphJson.g:3984:1: ( ']' )
            {
            // InternalElkGraphJson.g:3984:1: ( ']' )
            // InternalElkGraphJson.g:3985:2: ']'
            {
             before(grammarAccess.getElkEdgeTargetsAccess().getRightSquareBracketKeyword_3()); 
            match(input,63,FOLLOW_2); 
             after(grammarAccess.getElkEdgeTargetsAccess().getRightSquareBracketKeyword_3()); 

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
    // $ANTLR end "rule__ElkEdgeTargets__Group__3__Impl"


    // $ANTLR start "rule__ElkEdgeTargets__Group_1__0"
    // InternalElkGraphJson.g:3995:1: rule__ElkEdgeTargets__Group_1__0 : rule__ElkEdgeTargets__Group_1__0__Impl rule__ElkEdgeTargets__Group_1__1 ;
    public final void rule__ElkEdgeTargets__Group_1__0() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:3999:1: ( rule__ElkEdgeTargets__Group_1__0__Impl rule__ElkEdgeTargets__Group_1__1 )
            // InternalElkGraphJson.g:4000:2: rule__ElkEdgeTargets__Group_1__0__Impl rule__ElkEdgeTargets__Group_1__1
            {
            pushFollow(FOLLOW_5);
            rule__ElkEdgeTargets__Group_1__0__Impl();

            state._fsp--;

            pushFollow(FOLLOW_2);
            rule__ElkEdgeTargets__Group_1__1();

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
    // $ANTLR end "rule__ElkEdgeTargets__Group_1__0"


    // $ANTLR start "rule__ElkEdgeTargets__Group_1__0__Impl"
    // InternalElkGraphJson.g:4007:1: rule__ElkEdgeTargets__Group_1__0__Impl : ( ( rule__ElkEdgeTargets__TargetsAssignment_1_0 ) ) ;
    public final void rule__ElkEdgeTargets__Group_1__0__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:4011:1: ( ( ( rule__ElkEdgeTargets__TargetsAssignment_1_0 ) ) )
            // InternalElkGraphJson.g:4012:1: ( ( rule__ElkEdgeTargets__TargetsAssignment_1_0 ) )
            {
            // InternalElkGraphJson.g:4012:1: ( ( rule__ElkEdgeTargets__TargetsAssignment_1_0 ) )
            // InternalElkGraphJson.g:4013:2: ( rule__ElkEdgeTargets__TargetsAssignment_1_0 )
            {
             before(grammarAccess.getElkEdgeTargetsAccess().getTargetsAssignment_1_0()); 
            // InternalElkGraphJson.g:4014:2: ( rule__ElkEdgeTargets__TargetsAssignment_1_0 )
            // InternalElkGraphJson.g:4014:3: rule__ElkEdgeTargets__TargetsAssignment_1_0
            {
            pushFollow(FOLLOW_2);
            rule__ElkEdgeTargets__TargetsAssignment_1_0();

            state._fsp--;


            }

             after(grammarAccess.getElkEdgeTargetsAccess().getTargetsAssignment_1_0()); 

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
    // $ANTLR end "rule__ElkEdgeTargets__Group_1__0__Impl"


    // $ANTLR start "rule__ElkEdgeTargets__Group_1__1"
    // InternalElkGraphJson.g:4022:1: rule__ElkEdgeTargets__Group_1__1 : rule__ElkEdgeTargets__Group_1__1__Impl ;
    public final void rule__ElkEdgeTargets__Group_1__1() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:4026:1: ( rule__ElkEdgeTargets__Group_1__1__Impl )
            // InternalElkGraphJson.g:4027:2: rule__ElkEdgeTargets__Group_1__1__Impl
            {
            pushFollow(FOLLOW_2);
            rule__ElkEdgeTargets__Group_1__1__Impl();

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
    // $ANTLR end "rule__ElkEdgeTargets__Group_1__1"


    // $ANTLR start "rule__ElkEdgeTargets__Group_1__1__Impl"
    // InternalElkGraphJson.g:4033:1: rule__ElkEdgeTargets__Group_1__1__Impl : ( ( rule__ElkEdgeTargets__Group_1_1__0 )* ) ;
    public final void rule__ElkEdgeTargets__Group_1__1__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:4037:1: ( ( ( rule__ElkEdgeTargets__Group_1_1__0 )* ) )
            // InternalElkGraphJson.g:4038:1: ( ( rule__ElkEdgeTargets__Group_1_1__0 )* )
            {
            // InternalElkGraphJson.g:4038:1: ( ( rule__ElkEdgeTargets__Group_1_1__0 )* )
            // InternalElkGraphJson.g:4039:2: ( rule__ElkEdgeTargets__Group_1_1__0 )*
            {
             before(grammarAccess.getElkEdgeTargetsAccess().getGroup_1_1()); 
            // InternalElkGraphJson.g:4040:2: ( rule__ElkEdgeTargets__Group_1_1__0 )*
            loop40:
            do {
                int alt40=2;
                int LA40_0 = input.LA(1);

                if ( (LA40_0==59) ) {
                    int LA40_1 = input.LA(2);

                    if ( (LA40_1==RULE_STRING) ) {
                        alt40=1;
                    }


                }


                switch (alt40) {
            	case 1 :
            	    // InternalElkGraphJson.g:4040:3: rule__ElkEdgeTargets__Group_1_1__0
            	    {
            	    pushFollow(FOLLOW_6);
            	    rule__ElkEdgeTargets__Group_1_1__0();

            	    state._fsp--;


            	    }
            	    break;

            	default :
            	    break loop40;
                }
            } while (true);

             after(grammarAccess.getElkEdgeTargetsAccess().getGroup_1_1()); 

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
    // $ANTLR end "rule__ElkEdgeTargets__Group_1__1__Impl"


    // $ANTLR start "rule__ElkEdgeTargets__Group_1_1__0"
    // InternalElkGraphJson.g:4049:1: rule__ElkEdgeTargets__Group_1_1__0 : rule__ElkEdgeTargets__Group_1_1__0__Impl rule__ElkEdgeTargets__Group_1_1__1 ;
    public final void rule__ElkEdgeTargets__Group_1_1__0() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:4053:1: ( rule__ElkEdgeTargets__Group_1_1__0__Impl rule__ElkEdgeTargets__Group_1_1__1 )
            // InternalElkGraphJson.g:4054:2: rule__ElkEdgeTargets__Group_1_1__0__Impl rule__ElkEdgeTargets__Group_1_1__1
            {
            pushFollow(FOLLOW_12);
            rule__ElkEdgeTargets__Group_1_1__0__Impl();

            state._fsp--;

            pushFollow(FOLLOW_2);
            rule__ElkEdgeTargets__Group_1_1__1();

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
    // $ANTLR end "rule__ElkEdgeTargets__Group_1_1__0"


    // $ANTLR start "rule__ElkEdgeTargets__Group_1_1__0__Impl"
    // InternalElkGraphJson.g:4061:1: rule__ElkEdgeTargets__Group_1_1__0__Impl : ( ',' ) ;
    public final void rule__ElkEdgeTargets__Group_1_1__0__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:4065:1: ( ( ',' ) )
            // InternalElkGraphJson.g:4066:1: ( ',' )
            {
            // InternalElkGraphJson.g:4066:1: ( ',' )
            // InternalElkGraphJson.g:4067:2: ','
            {
             before(grammarAccess.getElkEdgeTargetsAccess().getCommaKeyword_1_1_0()); 
            match(input,59,FOLLOW_2); 
             after(grammarAccess.getElkEdgeTargetsAccess().getCommaKeyword_1_1_0()); 

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
    // $ANTLR end "rule__ElkEdgeTargets__Group_1_1__0__Impl"


    // $ANTLR start "rule__ElkEdgeTargets__Group_1_1__1"
    // InternalElkGraphJson.g:4076:1: rule__ElkEdgeTargets__Group_1_1__1 : rule__ElkEdgeTargets__Group_1_1__1__Impl ;
    public final void rule__ElkEdgeTargets__Group_1_1__1() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:4080:1: ( rule__ElkEdgeTargets__Group_1_1__1__Impl )
            // InternalElkGraphJson.g:4081:2: rule__ElkEdgeTargets__Group_1_1__1__Impl
            {
            pushFollow(FOLLOW_2);
            rule__ElkEdgeTargets__Group_1_1__1__Impl();

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
    // $ANTLR end "rule__ElkEdgeTargets__Group_1_1__1"


    // $ANTLR start "rule__ElkEdgeTargets__Group_1_1__1__Impl"
    // InternalElkGraphJson.g:4087:1: rule__ElkEdgeTargets__Group_1_1__1__Impl : ( ( rule__ElkEdgeTargets__TargetsAssignment_1_1_1 ) ) ;
    public final void rule__ElkEdgeTargets__Group_1_1__1__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:4091:1: ( ( ( rule__ElkEdgeTargets__TargetsAssignment_1_1_1 ) ) )
            // InternalElkGraphJson.g:4092:1: ( ( rule__ElkEdgeTargets__TargetsAssignment_1_1_1 ) )
            {
            // InternalElkGraphJson.g:4092:1: ( ( rule__ElkEdgeTargets__TargetsAssignment_1_1_1 ) )
            // InternalElkGraphJson.g:4093:2: ( rule__ElkEdgeTargets__TargetsAssignment_1_1_1 )
            {
             before(grammarAccess.getElkEdgeTargetsAccess().getTargetsAssignment_1_1_1()); 
            // InternalElkGraphJson.g:4094:2: ( rule__ElkEdgeTargets__TargetsAssignment_1_1_1 )
            // InternalElkGraphJson.g:4094:3: rule__ElkEdgeTargets__TargetsAssignment_1_1_1
            {
            pushFollow(FOLLOW_2);
            rule__ElkEdgeTargets__TargetsAssignment_1_1_1();

            state._fsp--;


            }

             after(grammarAccess.getElkEdgeTargetsAccess().getTargetsAssignment_1_1_1()); 

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
    // $ANTLR end "rule__ElkEdgeTargets__Group_1_1__1__Impl"


    // $ANTLR start "rule__ElkId__Group__0"
    // InternalElkGraphJson.g:4103:1: rule__ElkId__Group__0 : rule__ElkId__Group__0__Impl rule__ElkId__Group__1 ;
    public final void rule__ElkId__Group__0() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:4107:1: ( rule__ElkId__Group__0__Impl rule__ElkId__Group__1 )
            // InternalElkGraphJson.g:4108:2: rule__ElkId__Group__0__Impl rule__ElkId__Group__1
            {
            pushFollow(FOLLOW_8);
            rule__ElkId__Group__0__Impl();

            state._fsp--;

            pushFollow(FOLLOW_2);
            rule__ElkId__Group__1();

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
    // $ANTLR end "rule__ElkId__Group__0"


    // $ANTLR start "rule__ElkId__Group__0__Impl"
    // InternalElkGraphJson.g:4115:1: rule__ElkId__Group__0__Impl : ( ruleKeyId ) ;
    public final void rule__ElkId__Group__0__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:4119:1: ( ( ruleKeyId ) )
            // InternalElkGraphJson.g:4120:1: ( ruleKeyId )
            {
            // InternalElkGraphJson.g:4120:1: ( ruleKeyId )
            // InternalElkGraphJson.g:4121:2: ruleKeyId
            {
             before(grammarAccess.getElkIdAccess().getKeyIdParserRuleCall_0()); 
            pushFollow(FOLLOW_2);
            ruleKeyId();

            state._fsp--;

             after(grammarAccess.getElkIdAccess().getKeyIdParserRuleCall_0()); 

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
    // $ANTLR end "rule__ElkId__Group__0__Impl"


    // $ANTLR start "rule__ElkId__Group__1"
    // InternalElkGraphJson.g:4130:1: rule__ElkId__Group__1 : rule__ElkId__Group__1__Impl rule__ElkId__Group__2 ;
    public final void rule__ElkId__Group__1() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:4134:1: ( rule__ElkId__Group__1__Impl rule__ElkId__Group__2 )
            // InternalElkGraphJson.g:4135:2: rule__ElkId__Group__1__Impl rule__ElkId__Group__2
            {
            pushFollow(FOLLOW_12);
            rule__ElkId__Group__1__Impl();

            state._fsp--;

            pushFollow(FOLLOW_2);
            rule__ElkId__Group__2();

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
    // $ANTLR end "rule__ElkId__Group__1"


    // $ANTLR start "rule__ElkId__Group__1__Impl"
    // InternalElkGraphJson.g:4142:1: rule__ElkId__Group__1__Impl : ( ':' ) ;
    public final void rule__ElkId__Group__1__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:4146:1: ( ( ':' ) )
            // InternalElkGraphJson.g:4147:1: ( ':' )
            {
            // InternalElkGraphJson.g:4147:1: ( ':' )
            // InternalElkGraphJson.g:4148:2: ':'
            {
             before(grammarAccess.getElkIdAccess().getColonKeyword_1()); 
            match(input,61,FOLLOW_2); 
             after(grammarAccess.getElkIdAccess().getColonKeyword_1()); 

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
    // $ANTLR end "rule__ElkId__Group__1__Impl"


    // $ANTLR start "rule__ElkId__Group__2"
    // InternalElkGraphJson.g:4157:1: rule__ElkId__Group__2 : rule__ElkId__Group__2__Impl ;
    public final void rule__ElkId__Group__2() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:4161:1: ( rule__ElkId__Group__2__Impl )
            // InternalElkGraphJson.g:4162:2: rule__ElkId__Group__2__Impl
            {
            pushFollow(FOLLOW_2);
            rule__ElkId__Group__2__Impl();

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
    // $ANTLR end "rule__ElkId__Group__2"


    // $ANTLR start "rule__ElkId__Group__2__Impl"
    // InternalElkGraphJson.g:4168:1: rule__ElkId__Group__2__Impl : ( ( rule__ElkId__IdentifierAssignment_2 ) ) ;
    public final void rule__ElkId__Group__2__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:4172:1: ( ( ( rule__ElkId__IdentifierAssignment_2 ) ) )
            // InternalElkGraphJson.g:4173:1: ( ( rule__ElkId__IdentifierAssignment_2 ) )
            {
            // InternalElkGraphJson.g:4173:1: ( ( rule__ElkId__IdentifierAssignment_2 ) )
            // InternalElkGraphJson.g:4174:2: ( rule__ElkId__IdentifierAssignment_2 )
            {
             before(grammarAccess.getElkIdAccess().getIdentifierAssignment_2()); 
            // InternalElkGraphJson.g:4175:2: ( rule__ElkId__IdentifierAssignment_2 )
            // InternalElkGraphJson.g:4175:3: rule__ElkId__IdentifierAssignment_2
            {
            pushFollow(FOLLOW_2);
            rule__ElkId__IdentifierAssignment_2();

            state._fsp--;


            }

             after(grammarAccess.getElkIdAccess().getIdentifierAssignment_2()); 

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
    // $ANTLR end "rule__ElkId__Group__2__Impl"


    // $ANTLR start "rule__ElkNodeChildren__Group__0"
    // InternalElkGraphJson.g:4184:1: rule__ElkNodeChildren__Group__0 : rule__ElkNodeChildren__Group__0__Impl rule__ElkNodeChildren__Group__1 ;
    public final void rule__ElkNodeChildren__Group__0() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:4188:1: ( rule__ElkNodeChildren__Group__0__Impl rule__ElkNodeChildren__Group__1 )
            // InternalElkGraphJson.g:4189:2: rule__ElkNodeChildren__Group__0__Impl rule__ElkNodeChildren__Group__1
            {
            pushFollow(FOLLOW_15);
            rule__ElkNodeChildren__Group__0__Impl();

            state._fsp--;

            pushFollow(FOLLOW_2);
            rule__ElkNodeChildren__Group__1();

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
    // $ANTLR end "rule__ElkNodeChildren__Group__0"


    // $ANTLR start "rule__ElkNodeChildren__Group__0__Impl"
    // InternalElkGraphJson.g:4196:1: rule__ElkNodeChildren__Group__0__Impl : ( '[' ) ;
    public final void rule__ElkNodeChildren__Group__0__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:4200:1: ( ( '[' ) )
            // InternalElkGraphJson.g:4201:1: ( '[' )
            {
            // InternalElkGraphJson.g:4201:1: ( '[' )
            // InternalElkGraphJson.g:4202:2: '['
            {
             before(grammarAccess.getElkNodeChildrenAccess().getLeftSquareBracketKeyword_0()); 
            match(input,62,FOLLOW_2); 
             after(grammarAccess.getElkNodeChildrenAccess().getLeftSquareBracketKeyword_0()); 

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
    // $ANTLR end "rule__ElkNodeChildren__Group__0__Impl"


    // $ANTLR start "rule__ElkNodeChildren__Group__1"
    // InternalElkGraphJson.g:4211:1: rule__ElkNodeChildren__Group__1 : rule__ElkNodeChildren__Group__1__Impl rule__ElkNodeChildren__Group__2 ;
    public final void rule__ElkNodeChildren__Group__1() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:4215:1: ( rule__ElkNodeChildren__Group__1__Impl rule__ElkNodeChildren__Group__2 )
            // InternalElkGraphJson.g:4216:2: rule__ElkNodeChildren__Group__1__Impl rule__ElkNodeChildren__Group__2
            {
            pushFollow(FOLLOW_15);
            rule__ElkNodeChildren__Group__1__Impl();

            state._fsp--;

            pushFollow(FOLLOW_2);
            rule__ElkNodeChildren__Group__2();

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
    // $ANTLR end "rule__ElkNodeChildren__Group__1"


    // $ANTLR start "rule__ElkNodeChildren__Group__1__Impl"
    // InternalElkGraphJson.g:4223:1: rule__ElkNodeChildren__Group__1__Impl : ( ( rule__ElkNodeChildren__Group_1__0 )? ) ;
    public final void rule__ElkNodeChildren__Group__1__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:4227:1: ( ( ( rule__ElkNodeChildren__Group_1__0 )? ) )
            // InternalElkGraphJson.g:4228:1: ( ( rule__ElkNodeChildren__Group_1__0 )? )
            {
            // InternalElkGraphJson.g:4228:1: ( ( rule__ElkNodeChildren__Group_1__0 )? )
            // InternalElkGraphJson.g:4229:2: ( rule__ElkNodeChildren__Group_1__0 )?
            {
             before(grammarAccess.getElkNodeChildrenAccess().getGroup_1()); 
            // InternalElkGraphJson.g:4230:2: ( rule__ElkNodeChildren__Group_1__0 )?
            int alt41=2;
            int LA41_0 = input.LA(1);

            if ( (LA41_0==58) ) {
                alt41=1;
            }
            switch (alt41) {
                case 1 :
                    // InternalElkGraphJson.g:4230:3: rule__ElkNodeChildren__Group_1__0
                    {
                    pushFollow(FOLLOW_2);
                    rule__ElkNodeChildren__Group_1__0();

                    state._fsp--;


                    }
                    break;

            }

             after(grammarAccess.getElkNodeChildrenAccess().getGroup_1()); 

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
    // $ANTLR end "rule__ElkNodeChildren__Group__1__Impl"


    // $ANTLR start "rule__ElkNodeChildren__Group__2"
    // InternalElkGraphJson.g:4238:1: rule__ElkNodeChildren__Group__2 : rule__ElkNodeChildren__Group__2__Impl rule__ElkNodeChildren__Group__3 ;
    public final void rule__ElkNodeChildren__Group__2() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:4242:1: ( rule__ElkNodeChildren__Group__2__Impl rule__ElkNodeChildren__Group__3 )
            // InternalElkGraphJson.g:4243:2: rule__ElkNodeChildren__Group__2__Impl rule__ElkNodeChildren__Group__3
            {
            pushFollow(FOLLOW_15);
            rule__ElkNodeChildren__Group__2__Impl();

            state._fsp--;

            pushFollow(FOLLOW_2);
            rule__ElkNodeChildren__Group__3();

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
    // $ANTLR end "rule__ElkNodeChildren__Group__2"


    // $ANTLR start "rule__ElkNodeChildren__Group__2__Impl"
    // InternalElkGraphJson.g:4250:1: rule__ElkNodeChildren__Group__2__Impl : ( ( ',' )? ) ;
    public final void rule__ElkNodeChildren__Group__2__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:4254:1: ( ( ( ',' )? ) )
            // InternalElkGraphJson.g:4255:1: ( ( ',' )? )
            {
            // InternalElkGraphJson.g:4255:1: ( ( ',' )? )
            // InternalElkGraphJson.g:4256:2: ( ',' )?
            {
             before(grammarAccess.getElkNodeChildrenAccess().getCommaKeyword_2()); 
            // InternalElkGraphJson.g:4257:2: ( ',' )?
            int alt42=2;
            int LA42_0 = input.LA(1);

            if ( (LA42_0==59) ) {
                alt42=1;
            }
            switch (alt42) {
                case 1 :
                    // InternalElkGraphJson.g:4257:3: ','
                    {
                    match(input,59,FOLLOW_2); 

                    }
                    break;

            }

             after(grammarAccess.getElkNodeChildrenAccess().getCommaKeyword_2()); 

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
    // $ANTLR end "rule__ElkNodeChildren__Group__2__Impl"


    // $ANTLR start "rule__ElkNodeChildren__Group__3"
    // InternalElkGraphJson.g:4265:1: rule__ElkNodeChildren__Group__3 : rule__ElkNodeChildren__Group__3__Impl ;
    public final void rule__ElkNodeChildren__Group__3() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:4269:1: ( rule__ElkNodeChildren__Group__3__Impl )
            // InternalElkGraphJson.g:4270:2: rule__ElkNodeChildren__Group__3__Impl
            {
            pushFollow(FOLLOW_2);
            rule__ElkNodeChildren__Group__3__Impl();

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
    // $ANTLR end "rule__ElkNodeChildren__Group__3"


    // $ANTLR start "rule__ElkNodeChildren__Group__3__Impl"
    // InternalElkGraphJson.g:4276:1: rule__ElkNodeChildren__Group__3__Impl : ( ']' ) ;
    public final void rule__ElkNodeChildren__Group__3__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:4280:1: ( ( ']' ) )
            // InternalElkGraphJson.g:4281:1: ( ']' )
            {
            // InternalElkGraphJson.g:4281:1: ( ']' )
            // InternalElkGraphJson.g:4282:2: ']'
            {
             before(grammarAccess.getElkNodeChildrenAccess().getRightSquareBracketKeyword_3()); 
            match(input,63,FOLLOW_2); 
             after(grammarAccess.getElkNodeChildrenAccess().getRightSquareBracketKeyword_3()); 

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
    // $ANTLR end "rule__ElkNodeChildren__Group__3__Impl"


    // $ANTLR start "rule__ElkNodeChildren__Group_1__0"
    // InternalElkGraphJson.g:4292:1: rule__ElkNodeChildren__Group_1__0 : rule__ElkNodeChildren__Group_1__0__Impl rule__ElkNodeChildren__Group_1__1 ;
    public final void rule__ElkNodeChildren__Group_1__0() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:4296:1: ( rule__ElkNodeChildren__Group_1__0__Impl rule__ElkNodeChildren__Group_1__1 )
            // InternalElkGraphJson.g:4297:2: rule__ElkNodeChildren__Group_1__0__Impl rule__ElkNodeChildren__Group_1__1
            {
            pushFollow(FOLLOW_5);
            rule__ElkNodeChildren__Group_1__0__Impl();

            state._fsp--;

            pushFollow(FOLLOW_2);
            rule__ElkNodeChildren__Group_1__1();

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
    // $ANTLR end "rule__ElkNodeChildren__Group_1__0"


    // $ANTLR start "rule__ElkNodeChildren__Group_1__0__Impl"
    // InternalElkGraphJson.g:4304:1: rule__ElkNodeChildren__Group_1__0__Impl : ( ( rule__ElkNodeChildren__ChildrenAssignment_1_0 ) ) ;
    public final void rule__ElkNodeChildren__Group_1__0__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:4308:1: ( ( ( rule__ElkNodeChildren__ChildrenAssignment_1_0 ) ) )
            // InternalElkGraphJson.g:4309:1: ( ( rule__ElkNodeChildren__ChildrenAssignment_1_0 ) )
            {
            // InternalElkGraphJson.g:4309:1: ( ( rule__ElkNodeChildren__ChildrenAssignment_1_0 ) )
            // InternalElkGraphJson.g:4310:2: ( rule__ElkNodeChildren__ChildrenAssignment_1_0 )
            {
             before(grammarAccess.getElkNodeChildrenAccess().getChildrenAssignment_1_0()); 
            // InternalElkGraphJson.g:4311:2: ( rule__ElkNodeChildren__ChildrenAssignment_1_0 )
            // InternalElkGraphJson.g:4311:3: rule__ElkNodeChildren__ChildrenAssignment_1_0
            {
            pushFollow(FOLLOW_2);
            rule__ElkNodeChildren__ChildrenAssignment_1_0();

            state._fsp--;


            }

             after(grammarAccess.getElkNodeChildrenAccess().getChildrenAssignment_1_0()); 

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
    // $ANTLR end "rule__ElkNodeChildren__Group_1__0__Impl"


    // $ANTLR start "rule__ElkNodeChildren__Group_1__1"
    // InternalElkGraphJson.g:4319:1: rule__ElkNodeChildren__Group_1__1 : rule__ElkNodeChildren__Group_1__1__Impl ;
    public final void rule__ElkNodeChildren__Group_1__1() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:4323:1: ( rule__ElkNodeChildren__Group_1__1__Impl )
            // InternalElkGraphJson.g:4324:2: rule__ElkNodeChildren__Group_1__1__Impl
            {
            pushFollow(FOLLOW_2);
            rule__ElkNodeChildren__Group_1__1__Impl();

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
    // $ANTLR end "rule__ElkNodeChildren__Group_1__1"


    // $ANTLR start "rule__ElkNodeChildren__Group_1__1__Impl"
    // InternalElkGraphJson.g:4330:1: rule__ElkNodeChildren__Group_1__1__Impl : ( ( rule__ElkNodeChildren__Group_1_1__0 )* ) ;
    public final void rule__ElkNodeChildren__Group_1__1__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:4334:1: ( ( ( rule__ElkNodeChildren__Group_1_1__0 )* ) )
            // InternalElkGraphJson.g:4335:1: ( ( rule__ElkNodeChildren__Group_1_1__0 )* )
            {
            // InternalElkGraphJson.g:4335:1: ( ( rule__ElkNodeChildren__Group_1_1__0 )* )
            // InternalElkGraphJson.g:4336:2: ( rule__ElkNodeChildren__Group_1_1__0 )*
            {
             before(grammarAccess.getElkNodeChildrenAccess().getGroup_1_1()); 
            // InternalElkGraphJson.g:4337:2: ( rule__ElkNodeChildren__Group_1_1__0 )*
            loop43:
            do {
                int alt43=2;
                int LA43_0 = input.LA(1);

                if ( (LA43_0==59) ) {
                    int LA43_1 = input.LA(2);

                    if ( (LA43_1==58) ) {
                        alt43=1;
                    }


                }


                switch (alt43) {
            	case 1 :
            	    // InternalElkGraphJson.g:4337:3: rule__ElkNodeChildren__Group_1_1__0
            	    {
            	    pushFollow(FOLLOW_6);
            	    rule__ElkNodeChildren__Group_1_1__0();

            	    state._fsp--;


            	    }
            	    break;

            	default :
            	    break loop43;
                }
            } while (true);

             after(grammarAccess.getElkNodeChildrenAccess().getGroup_1_1()); 

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
    // $ANTLR end "rule__ElkNodeChildren__Group_1__1__Impl"


    // $ANTLR start "rule__ElkNodeChildren__Group_1_1__0"
    // InternalElkGraphJson.g:4346:1: rule__ElkNodeChildren__Group_1_1__0 : rule__ElkNodeChildren__Group_1_1__0__Impl rule__ElkNodeChildren__Group_1_1__1 ;
    public final void rule__ElkNodeChildren__Group_1_1__0() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:4350:1: ( rule__ElkNodeChildren__Group_1_1__0__Impl rule__ElkNodeChildren__Group_1_1__1 )
            // InternalElkGraphJson.g:4351:2: rule__ElkNodeChildren__Group_1_1__0__Impl rule__ElkNodeChildren__Group_1_1__1
            {
            pushFollow(FOLLOW_3);
            rule__ElkNodeChildren__Group_1_1__0__Impl();

            state._fsp--;

            pushFollow(FOLLOW_2);
            rule__ElkNodeChildren__Group_1_1__1();

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
    // $ANTLR end "rule__ElkNodeChildren__Group_1_1__0"


    // $ANTLR start "rule__ElkNodeChildren__Group_1_1__0__Impl"
    // InternalElkGraphJson.g:4358:1: rule__ElkNodeChildren__Group_1_1__0__Impl : ( ',' ) ;
    public final void rule__ElkNodeChildren__Group_1_1__0__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:4362:1: ( ( ',' ) )
            // InternalElkGraphJson.g:4363:1: ( ',' )
            {
            // InternalElkGraphJson.g:4363:1: ( ',' )
            // InternalElkGraphJson.g:4364:2: ','
            {
             before(grammarAccess.getElkNodeChildrenAccess().getCommaKeyword_1_1_0()); 
            match(input,59,FOLLOW_2); 
             after(grammarAccess.getElkNodeChildrenAccess().getCommaKeyword_1_1_0()); 

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
    // $ANTLR end "rule__ElkNodeChildren__Group_1_1__0__Impl"


    // $ANTLR start "rule__ElkNodeChildren__Group_1_1__1"
    // InternalElkGraphJson.g:4373:1: rule__ElkNodeChildren__Group_1_1__1 : rule__ElkNodeChildren__Group_1_1__1__Impl ;
    public final void rule__ElkNodeChildren__Group_1_1__1() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:4377:1: ( rule__ElkNodeChildren__Group_1_1__1__Impl )
            // InternalElkGraphJson.g:4378:2: rule__ElkNodeChildren__Group_1_1__1__Impl
            {
            pushFollow(FOLLOW_2);
            rule__ElkNodeChildren__Group_1_1__1__Impl();

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
    // $ANTLR end "rule__ElkNodeChildren__Group_1_1__1"


    // $ANTLR start "rule__ElkNodeChildren__Group_1_1__1__Impl"
    // InternalElkGraphJson.g:4384:1: rule__ElkNodeChildren__Group_1_1__1__Impl : ( ( rule__ElkNodeChildren__ChildrenAssignment_1_1_1 ) ) ;
    public final void rule__ElkNodeChildren__Group_1_1__1__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:4388:1: ( ( ( rule__ElkNodeChildren__ChildrenAssignment_1_1_1 ) ) )
            // InternalElkGraphJson.g:4389:1: ( ( rule__ElkNodeChildren__ChildrenAssignment_1_1_1 ) )
            {
            // InternalElkGraphJson.g:4389:1: ( ( rule__ElkNodeChildren__ChildrenAssignment_1_1_1 ) )
            // InternalElkGraphJson.g:4390:2: ( rule__ElkNodeChildren__ChildrenAssignment_1_1_1 )
            {
             before(grammarAccess.getElkNodeChildrenAccess().getChildrenAssignment_1_1_1()); 
            // InternalElkGraphJson.g:4391:2: ( rule__ElkNodeChildren__ChildrenAssignment_1_1_1 )
            // InternalElkGraphJson.g:4391:3: rule__ElkNodeChildren__ChildrenAssignment_1_1_1
            {
            pushFollow(FOLLOW_2);
            rule__ElkNodeChildren__ChildrenAssignment_1_1_1();

            state._fsp--;


            }

             after(grammarAccess.getElkNodeChildrenAccess().getChildrenAssignment_1_1_1()); 

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
    // $ANTLR end "rule__ElkNodeChildren__Group_1_1__1__Impl"


    // $ANTLR start "rule__ElkNodePorts__Group__0"
    // InternalElkGraphJson.g:4400:1: rule__ElkNodePorts__Group__0 : rule__ElkNodePorts__Group__0__Impl rule__ElkNodePorts__Group__1 ;
    public final void rule__ElkNodePorts__Group__0() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:4404:1: ( rule__ElkNodePorts__Group__0__Impl rule__ElkNodePorts__Group__1 )
            // InternalElkGraphJson.g:4405:2: rule__ElkNodePorts__Group__0__Impl rule__ElkNodePorts__Group__1
            {
            pushFollow(FOLLOW_15);
            rule__ElkNodePorts__Group__0__Impl();

            state._fsp--;

            pushFollow(FOLLOW_2);
            rule__ElkNodePorts__Group__1();

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
    // $ANTLR end "rule__ElkNodePorts__Group__0"


    // $ANTLR start "rule__ElkNodePorts__Group__0__Impl"
    // InternalElkGraphJson.g:4412:1: rule__ElkNodePorts__Group__0__Impl : ( '[' ) ;
    public final void rule__ElkNodePorts__Group__0__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:4416:1: ( ( '[' ) )
            // InternalElkGraphJson.g:4417:1: ( '[' )
            {
            // InternalElkGraphJson.g:4417:1: ( '[' )
            // InternalElkGraphJson.g:4418:2: '['
            {
             before(grammarAccess.getElkNodePortsAccess().getLeftSquareBracketKeyword_0()); 
            match(input,62,FOLLOW_2); 
             after(grammarAccess.getElkNodePortsAccess().getLeftSquareBracketKeyword_0()); 

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
    // $ANTLR end "rule__ElkNodePorts__Group__0__Impl"


    // $ANTLR start "rule__ElkNodePorts__Group__1"
    // InternalElkGraphJson.g:4427:1: rule__ElkNodePorts__Group__1 : rule__ElkNodePorts__Group__1__Impl rule__ElkNodePorts__Group__2 ;
    public final void rule__ElkNodePorts__Group__1() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:4431:1: ( rule__ElkNodePorts__Group__1__Impl rule__ElkNodePorts__Group__2 )
            // InternalElkGraphJson.g:4432:2: rule__ElkNodePorts__Group__1__Impl rule__ElkNodePorts__Group__2
            {
            pushFollow(FOLLOW_15);
            rule__ElkNodePorts__Group__1__Impl();

            state._fsp--;

            pushFollow(FOLLOW_2);
            rule__ElkNodePorts__Group__2();

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
    // $ANTLR end "rule__ElkNodePorts__Group__1"


    // $ANTLR start "rule__ElkNodePorts__Group__1__Impl"
    // InternalElkGraphJson.g:4439:1: rule__ElkNodePorts__Group__1__Impl : ( ( rule__ElkNodePorts__Group_1__0 )? ) ;
    public final void rule__ElkNodePorts__Group__1__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:4443:1: ( ( ( rule__ElkNodePorts__Group_1__0 )? ) )
            // InternalElkGraphJson.g:4444:1: ( ( rule__ElkNodePorts__Group_1__0 )? )
            {
            // InternalElkGraphJson.g:4444:1: ( ( rule__ElkNodePorts__Group_1__0 )? )
            // InternalElkGraphJson.g:4445:2: ( rule__ElkNodePorts__Group_1__0 )?
            {
             before(grammarAccess.getElkNodePortsAccess().getGroup_1()); 
            // InternalElkGraphJson.g:4446:2: ( rule__ElkNodePorts__Group_1__0 )?
            int alt44=2;
            int LA44_0 = input.LA(1);

            if ( (LA44_0==58) ) {
                alt44=1;
            }
            switch (alt44) {
                case 1 :
                    // InternalElkGraphJson.g:4446:3: rule__ElkNodePorts__Group_1__0
                    {
                    pushFollow(FOLLOW_2);
                    rule__ElkNodePorts__Group_1__0();

                    state._fsp--;


                    }
                    break;

            }

             after(grammarAccess.getElkNodePortsAccess().getGroup_1()); 

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
    // $ANTLR end "rule__ElkNodePorts__Group__1__Impl"


    // $ANTLR start "rule__ElkNodePorts__Group__2"
    // InternalElkGraphJson.g:4454:1: rule__ElkNodePorts__Group__2 : rule__ElkNodePorts__Group__2__Impl rule__ElkNodePorts__Group__3 ;
    public final void rule__ElkNodePorts__Group__2() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:4458:1: ( rule__ElkNodePorts__Group__2__Impl rule__ElkNodePorts__Group__3 )
            // InternalElkGraphJson.g:4459:2: rule__ElkNodePorts__Group__2__Impl rule__ElkNodePorts__Group__3
            {
            pushFollow(FOLLOW_15);
            rule__ElkNodePorts__Group__2__Impl();

            state._fsp--;

            pushFollow(FOLLOW_2);
            rule__ElkNodePorts__Group__3();

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
    // $ANTLR end "rule__ElkNodePorts__Group__2"


    // $ANTLR start "rule__ElkNodePorts__Group__2__Impl"
    // InternalElkGraphJson.g:4466:1: rule__ElkNodePorts__Group__2__Impl : ( ( ',' )? ) ;
    public final void rule__ElkNodePorts__Group__2__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:4470:1: ( ( ( ',' )? ) )
            // InternalElkGraphJson.g:4471:1: ( ( ',' )? )
            {
            // InternalElkGraphJson.g:4471:1: ( ( ',' )? )
            // InternalElkGraphJson.g:4472:2: ( ',' )?
            {
             before(grammarAccess.getElkNodePortsAccess().getCommaKeyword_2()); 
            // InternalElkGraphJson.g:4473:2: ( ',' )?
            int alt45=2;
            int LA45_0 = input.LA(1);

            if ( (LA45_0==59) ) {
                alt45=1;
            }
            switch (alt45) {
                case 1 :
                    // InternalElkGraphJson.g:4473:3: ','
                    {
                    match(input,59,FOLLOW_2); 

                    }
                    break;

            }

             after(grammarAccess.getElkNodePortsAccess().getCommaKeyword_2()); 

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
    // $ANTLR end "rule__ElkNodePorts__Group__2__Impl"


    // $ANTLR start "rule__ElkNodePorts__Group__3"
    // InternalElkGraphJson.g:4481:1: rule__ElkNodePorts__Group__3 : rule__ElkNodePorts__Group__3__Impl ;
    public final void rule__ElkNodePorts__Group__3() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:4485:1: ( rule__ElkNodePorts__Group__3__Impl )
            // InternalElkGraphJson.g:4486:2: rule__ElkNodePorts__Group__3__Impl
            {
            pushFollow(FOLLOW_2);
            rule__ElkNodePorts__Group__3__Impl();

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
    // $ANTLR end "rule__ElkNodePorts__Group__3"


    // $ANTLR start "rule__ElkNodePorts__Group__3__Impl"
    // InternalElkGraphJson.g:4492:1: rule__ElkNodePorts__Group__3__Impl : ( ']' ) ;
    public final void rule__ElkNodePorts__Group__3__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:4496:1: ( ( ']' ) )
            // InternalElkGraphJson.g:4497:1: ( ']' )
            {
            // InternalElkGraphJson.g:4497:1: ( ']' )
            // InternalElkGraphJson.g:4498:2: ']'
            {
             before(grammarAccess.getElkNodePortsAccess().getRightSquareBracketKeyword_3()); 
            match(input,63,FOLLOW_2); 
             after(grammarAccess.getElkNodePortsAccess().getRightSquareBracketKeyword_3()); 

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
    // $ANTLR end "rule__ElkNodePorts__Group__3__Impl"


    // $ANTLR start "rule__ElkNodePorts__Group_1__0"
    // InternalElkGraphJson.g:4508:1: rule__ElkNodePorts__Group_1__0 : rule__ElkNodePorts__Group_1__0__Impl rule__ElkNodePorts__Group_1__1 ;
    public final void rule__ElkNodePorts__Group_1__0() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:4512:1: ( rule__ElkNodePorts__Group_1__0__Impl rule__ElkNodePorts__Group_1__1 )
            // InternalElkGraphJson.g:4513:2: rule__ElkNodePorts__Group_1__0__Impl rule__ElkNodePorts__Group_1__1
            {
            pushFollow(FOLLOW_5);
            rule__ElkNodePorts__Group_1__0__Impl();

            state._fsp--;

            pushFollow(FOLLOW_2);
            rule__ElkNodePorts__Group_1__1();

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
    // $ANTLR end "rule__ElkNodePorts__Group_1__0"


    // $ANTLR start "rule__ElkNodePorts__Group_1__0__Impl"
    // InternalElkGraphJson.g:4520:1: rule__ElkNodePorts__Group_1__0__Impl : ( ( rule__ElkNodePorts__PortsAssignment_1_0 ) ) ;
    public final void rule__ElkNodePorts__Group_1__0__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:4524:1: ( ( ( rule__ElkNodePorts__PortsAssignment_1_0 ) ) )
            // InternalElkGraphJson.g:4525:1: ( ( rule__ElkNodePorts__PortsAssignment_1_0 ) )
            {
            // InternalElkGraphJson.g:4525:1: ( ( rule__ElkNodePorts__PortsAssignment_1_0 ) )
            // InternalElkGraphJson.g:4526:2: ( rule__ElkNodePorts__PortsAssignment_1_0 )
            {
             before(grammarAccess.getElkNodePortsAccess().getPortsAssignment_1_0()); 
            // InternalElkGraphJson.g:4527:2: ( rule__ElkNodePorts__PortsAssignment_1_0 )
            // InternalElkGraphJson.g:4527:3: rule__ElkNodePorts__PortsAssignment_1_0
            {
            pushFollow(FOLLOW_2);
            rule__ElkNodePorts__PortsAssignment_1_0();

            state._fsp--;


            }

             after(grammarAccess.getElkNodePortsAccess().getPortsAssignment_1_0()); 

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
    // $ANTLR end "rule__ElkNodePorts__Group_1__0__Impl"


    // $ANTLR start "rule__ElkNodePorts__Group_1__1"
    // InternalElkGraphJson.g:4535:1: rule__ElkNodePorts__Group_1__1 : rule__ElkNodePorts__Group_1__1__Impl ;
    public final void rule__ElkNodePorts__Group_1__1() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:4539:1: ( rule__ElkNodePorts__Group_1__1__Impl )
            // InternalElkGraphJson.g:4540:2: rule__ElkNodePorts__Group_1__1__Impl
            {
            pushFollow(FOLLOW_2);
            rule__ElkNodePorts__Group_1__1__Impl();

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
    // $ANTLR end "rule__ElkNodePorts__Group_1__1"


    // $ANTLR start "rule__ElkNodePorts__Group_1__1__Impl"
    // InternalElkGraphJson.g:4546:1: rule__ElkNodePorts__Group_1__1__Impl : ( ( rule__ElkNodePorts__Group_1_1__0 )* ) ;
    public final void rule__ElkNodePorts__Group_1__1__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:4550:1: ( ( ( rule__ElkNodePorts__Group_1_1__0 )* ) )
            // InternalElkGraphJson.g:4551:1: ( ( rule__ElkNodePorts__Group_1_1__0 )* )
            {
            // InternalElkGraphJson.g:4551:1: ( ( rule__ElkNodePorts__Group_1_1__0 )* )
            // InternalElkGraphJson.g:4552:2: ( rule__ElkNodePorts__Group_1_1__0 )*
            {
             before(grammarAccess.getElkNodePortsAccess().getGroup_1_1()); 
            // InternalElkGraphJson.g:4553:2: ( rule__ElkNodePorts__Group_1_1__0 )*
            loop46:
            do {
                int alt46=2;
                int LA46_0 = input.LA(1);

                if ( (LA46_0==59) ) {
                    int LA46_1 = input.LA(2);

                    if ( (LA46_1==58) ) {
                        alt46=1;
                    }


                }


                switch (alt46) {
            	case 1 :
            	    // InternalElkGraphJson.g:4553:3: rule__ElkNodePorts__Group_1_1__0
            	    {
            	    pushFollow(FOLLOW_6);
            	    rule__ElkNodePorts__Group_1_1__0();

            	    state._fsp--;


            	    }
            	    break;

            	default :
            	    break loop46;
                }
            } while (true);

             after(grammarAccess.getElkNodePortsAccess().getGroup_1_1()); 

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
    // $ANTLR end "rule__ElkNodePorts__Group_1__1__Impl"


    // $ANTLR start "rule__ElkNodePorts__Group_1_1__0"
    // InternalElkGraphJson.g:4562:1: rule__ElkNodePorts__Group_1_1__0 : rule__ElkNodePorts__Group_1_1__0__Impl rule__ElkNodePorts__Group_1_1__1 ;
    public final void rule__ElkNodePorts__Group_1_1__0() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:4566:1: ( rule__ElkNodePorts__Group_1_1__0__Impl rule__ElkNodePorts__Group_1_1__1 )
            // InternalElkGraphJson.g:4567:2: rule__ElkNodePorts__Group_1_1__0__Impl rule__ElkNodePorts__Group_1_1__1
            {
            pushFollow(FOLLOW_3);
            rule__ElkNodePorts__Group_1_1__0__Impl();

            state._fsp--;

            pushFollow(FOLLOW_2);
            rule__ElkNodePorts__Group_1_1__1();

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
    // $ANTLR end "rule__ElkNodePorts__Group_1_1__0"


    // $ANTLR start "rule__ElkNodePorts__Group_1_1__0__Impl"
    // InternalElkGraphJson.g:4574:1: rule__ElkNodePorts__Group_1_1__0__Impl : ( ',' ) ;
    public final void rule__ElkNodePorts__Group_1_1__0__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:4578:1: ( ( ',' ) )
            // InternalElkGraphJson.g:4579:1: ( ',' )
            {
            // InternalElkGraphJson.g:4579:1: ( ',' )
            // InternalElkGraphJson.g:4580:2: ','
            {
             before(grammarAccess.getElkNodePortsAccess().getCommaKeyword_1_1_0()); 
            match(input,59,FOLLOW_2); 
             after(grammarAccess.getElkNodePortsAccess().getCommaKeyword_1_1_0()); 

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
    // $ANTLR end "rule__ElkNodePorts__Group_1_1__0__Impl"


    // $ANTLR start "rule__ElkNodePorts__Group_1_1__1"
    // InternalElkGraphJson.g:4589:1: rule__ElkNodePorts__Group_1_1__1 : rule__ElkNodePorts__Group_1_1__1__Impl ;
    public final void rule__ElkNodePorts__Group_1_1__1() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:4593:1: ( rule__ElkNodePorts__Group_1_1__1__Impl )
            // InternalElkGraphJson.g:4594:2: rule__ElkNodePorts__Group_1_1__1__Impl
            {
            pushFollow(FOLLOW_2);
            rule__ElkNodePorts__Group_1_1__1__Impl();

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
    // $ANTLR end "rule__ElkNodePorts__Group_1_1__1"


    // $ANTLR start "rule__ElkNodePorts__Group_1_1__1__Impl"
    // InternalElkGraphJson.g:4600:1: rule__ElkNodePorts__Group_1_1__1__Impl : ( ( rule__ElkNodePorts__PortsAssignment_1_1_1 ) ) ;
    public final void rule__ElkNodePorts__Group_1_1__1__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:4604:1: ( ( ( rule__ElkNodePorts__PortsAssignment_1_1_1 ) ) )
            // InternalElkGraphJson.g:4605:1: ( ( rule__ElkNodePorts__PortsAssignment_1_1_1 ) )
            {
            // InternalElkGraphJson.g:4605:1: ( ( rule__ElkNodePorts__PortsAssignment_1_1_1 ) )
            // InternalElkGraphJson.g:4606:2: ( rule__ElkNodePorts__PortsAssignment_1_1_1 )
            {
             before(grammarAccess.getElkNodePortsAccess().getPortsAssignment_1_1_1()); 
            // InternalElkGraphJson.g:4607:2: ( rule__ElkNodePorts__PortsAssignment_1_1_1 )
            // InternalElkGraphJson.g:4607:3: rule__ElkNodePorts__PortsAssignment_1_1_1
            {
            pushFollow(FOLLOW_2);
            rule__ElkNodePorts__PortsAssignment_1_1_1();

            state._fsp--;


            }

             after(grammarAccess.getElkNodePortsAccess().getPortsAssignment_1_1_1()); 

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
    // $ANTLR end "rule__ElkNodePorts__Group_1_1__1__Impl"


    // $ANTLR start "rule__ElkNodeEdges__Group__0"
    // InternalElkGraphJson.g:4616:1: rule__ElkNodeEdges__Group__0 : rule__ElkNodeEdges__Group__0__Impl rule__ElkNodeEdges__Group__1 ;
    public final void rule__ElkNodeEdges__Group__0() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:4620:1: ( rule__ElkNodeEdges__Group__0__Impl rule__ElkNodeEdges__Group__1 )
            // InternalElkGraphJson.g:4621:2: rule__ElkNodeEdges__Group__0__Impl rule__ElkNodeEdges__Group__1
            {
            pushFollow(FOLLOW_15);
            rule__ElkNodeEdges__Group__0__Impl();

            state._fsp--;

            pushFollow(FOLLOW_2);
            rule__ElkNodeEdges__Group__1();

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
    // $ANTLR end "rule__ElkNodeEdges__Group__0"


    // $ANTLR start "rule__ElkNodeEdges__Group__0__Impl"
    // InternalElkGraphJson.g:4628:1: rule__ElkNodeEdges__Group__0__Impl : ( '[' ) ;
    public final void rule__ElkNodeEdges__Group__0__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:4632:1: ( ( '[' ) )
            // InternalElkGraphJson.g:4633:1: ( '[' )
            {
            // InternalElkGraphJson.g:4633:1: ( '[' )
            // InternalElkGraphJson.g:4634:2: '['
            {
             before(grammarAccess.getElkNodeEdgesAccess().getLeftSquareBracketKeyword_0()); 
            match(input,62,FOLLOW_2); 
             after(grammarAccess.getElkNodeEdgesAccess().getLeftSquareBracketKeyword_0()); 

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
    // $ANTLR end "rule__ElkNodeEdges__Group__0__Impl"


    // $ANTLR start "rule__ElkNodeEdges__Group__1"
    // InternalElkGraphJson.g:4643:1: rule__ElkNodeEdges__Group__1 : rule__ElkNodeEdges__Group__1__Impl rule__ElkNodeEdges__Group__2 ;
    public final void rule__ElkNodeEdges__Group__1() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:4647:1: ( rule__ElkNodeEdges__Group__1__Impl rule__ElkNodeEdges__Group__2 )
            // InternalElkGraphJson.g:4648:2: rule__ElkNodeEdges__Group__1__Impl rule__ElkNodeEdges__Group__2
            {
            pushFollow(FOLLOW_15);
            rule__ElkNodeEdges__Group__1__Impl();

            state._fsp--;

            pushFollow(FOLLOW_2);
            rule__ElkNodeEdges__Group__2();

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
    // $ANTLR end "rule__ElkNodeEdges__Group__1"


    // $ANTLR start "rule__ElkNodeEdges__Group__1__Impl"
    // InternalElkGraphJson.g:4655:1: rule__ElkNodeEdges__Group__1__Impl : ( ( rule__ElkNodeEdges__Group_1__0 )? ) ;
    public final void rule__ElkNodeEdges__Group__1__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:4659:1: ( ( ( rule__ElkNodeEdges__Group_1__0 )? ) )
            // InternalElkGraphJson.g:4660:1: ( ( rule__ElkNodeEdges__Group_1__0 )? )
            {
            // InternalElkGraphJson.g:4660:1: ( ( rule__ElkNodeEdges__Group_1__0 )? )
            // InternalElkGraphJson.g:4661:2: ( rule__ElkNodeEdges__Group_1__0 )?
            {
             before(grammarAccess.getElkNodeEdgesAccess().getGroup_1()); 
            // InternalElkGraphJson.g:4662:2: ( rule__ElkNodeEdges__Group_1__0 )?
            int alt47=2;
            int LA47_0 = input.LA(1);

            if ( (LA47_0==58) ) {
                alt47=1;
            }
            switch (alt47) {
                case 1 :
                    // InternalElkGraphJson.g:4662:3: rule__ElkNodeEdges__Group_1__0
                    {
                    pushFollow(FOLLOW_2);
                    rule__ElkNodeEdges__Group_1__0();

                    state._fsp--;


                    }
                    break;

            }

             after(grammarAccess.getElkNodeEdgesAccess().getGroup_1()); 

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
    // $ANTLR end "rule__ElkNodeEdges__Group__1__Impl"


    // $ANTLR start "rule__ElkNodeEdges__Group__2"
    // InternalElkGraphJson.g:4670:1: rule__ElkNodeEdges__Group__2 : rule__ElkNodeEdges__Group__2__Impl rule__ElkNodeEdges__Group__3 ;
    public final void rule__ElkNodeEdges__Group__2() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:4674:1: ( rule__ElkNodeEdges__Group__2__Impl rule__ElkNodeEdges__Group__3 )
            // InternalElkGraphJson.g:4675:2: rule__ElkNodeEdges__Group__2__Impl rule__ElkNodeEdges__Group__3
            {
            pushFollow(FOLLOW_15);
            rule__ElkNodeEdges__Group__2__Impl();

            state._fsp--;

            pushFollow(FOLLOW_2);
            rule__ElkNodeEdges__Group__3();

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
    // $ANTLR end "rule__ElkNodeEdges__Group__2"


    // $ANTLR start "rule__ElkNodeEdges__Group__2__Impl"
    // InternalElkGraphJson.g:4682:1: rule__ElkNodeEdges__Group__2__Impl : ( ( ',' )? ) ;
    public final void rule__ElkNodeEdges__Group__2__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:4686:1: ( ( ( ',' )? ) )
            // InternalElkGraphJson.g:4687:1: ( ( ',' )? )
            {
            // InternalElkGraphJson.g:4687:1: ( ( ',' )? )
            // InternalElkGraphJson.g:4688:2: ( ',' )?
            {
             before(grammarAccess.getElkNodeEdgesAccess().getCommaKeyword_2()); 
            // InternalElkGraphJson.g:4689:2: ( ',' )?
            int alt48=2;
            int LA48_0 = input.LA(1);

            if ( (LA48_0==59) ) {
                alt48=1;
            }
            switch (alt48) {
                case 1 :
                    // InternalElkGraphJson.g:4689:3: ','
                    {
                    match(input,59,FOLLOW_2); 

                    }
                    break;

            }

             after(grammarAccess.getElkNodeEdgesAccess().getCommaKeyword_2()); 

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
    // $ANTLR end "rule__ElkNodeEdges__Group__2__Impl"


    // $ANTLR start "rule__ElkNodeEdges__Group__3"
    // InternalElkGraphJson.g:4697:1: rule__ElkNodeEdges__Group__3 : rule__ElkNodeEdges__Group__3__Impl ;
    public final void rule__ElkNodeEdges__Group__3() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:4701:1: ( rule__ElkNodeEdges__Group__3__Impl )
            // InternalElkGraphJson.g:4702:2: rule__ElkNodeEdges__Group__3__Impl
            {
            pushFollow(FOLLOW_2);
            rule__ElkNodeEdges__Group__3__Impl();

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
    // $ANTLR end "rule__ElkNodeEdges__Group__3"


    // $ANTLR start "rule__ElkNodeEdges__Group__3__Impl"
    // InternalElkGraphJson.g:4708:1: rule__ElkNodeEdges__Group__3__Impl : ( ']' ) ;
    public final void rule__ElkNodeEdges__Group__3__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:4712:1: ( ( ']' ) )
            // InternalElkGraphJson.g:4713:1: ( ']' )
            {
            // InternalElkGraphJson.g:4713:1: ( ']' )
            // InternalElkGraphJson.g:4714:2: ']'
            {
             before(grammarAccess.getElkNodeEdgesAccess().getRightSquareBracketKeyword_3()); 
            match(input,63,FOLLOW_2); 
             after(grammarAccess.getElkNodeEdgesAccess().getRightSquareBracketKeyword_3()); 

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
    // $ANTLR end "rule__ElkNodeEdges__Group__3__Impl"


    // $ANTLR start "rule__ElkNodeEdges__Group_1__0"
    // InternalElkGraphJson.g:4724:1: rule__ElkNodeEdges__Group_1__0 : rule__ElkNodeEdges__Group_1__0__Impl rule__ElkNodeEdges__Group_1__1 ;
    public final void rule__ElkNodeEdges__Group_1__0() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:4728:1: ( rule__ElkNodeEdges__Group_1__0__Impl rule__ElkNodeEdges__Group_1__1 )
            // InternalElkGraphJson.g:4729:2: rule__ElkNodeEdges__Group_1__0__Impl rule__ElkNodeEdges__Group_1__1
            {
            pushFollow(FOLLOW_5);
            rule__ElkNodeEdges__Group_1__0__Impl();

            state._fsp--;

            pushFollow(FOLLOW_2);
            rule__ElkNodeEdges__Group_1__1();

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
    // $ANTLR end "rule__ElkNodeEdges__Group_1__0"


    // $ANTLR start "rule__ElkNodeEdges__Group_1__0__Impl"
    // InternalElkGraphJson.g:4736:1: rule__ElkNodeEdges__Group_1__0__Impl : ( ( rule__ElkNodeEdges__ContainedEdgesAssignment_1_0 ) ) ;
    public final void rule__ElkNodeEdges__Group_1__0__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:4740:1: ( ( ( rule__ElkNodeEdges__ContainedEdgesAssignment_1_0 ) ) )
            // InternalElkGraphJson.g:4741:1: ( ( rule__ElkNodeEdges__ContainedEdgesAssignment_1_0 ) )
            {
            // InternalElkGraphJson.g:4741:1: ( ( rule__ElkNodeEdges__ContainedEdgesAssignment_1_0 ) )
            // InternalElkGraphJson.g:4742:2: ( rule__ElkNodeEdges__ContainedEdgesAssignment_1_0 )
            {
             before(grammarAccess.getElkNodeEdgesAccess().getContainedEdgesAssignment_1_0()); 
            // InternalElkGraphJson.g:4743:2: ( rule__ElkNodeEdges__ContainedEdgesAssignment_1_0 )
            // InternalElkGraphJson.g:4743:3: rule__ElkNodeEdges__ContainedEdgesAssignment_1_0
            {
            pushFollow(FOLLOW_2);
            rule__ElkNodeEdges__ContainedEdgesAssignment_1_0();

            state._fsp--;


            }

             after(grammarAccess.getElkNodeEdgesAccess().getContainedEdgesAssignment_1_0()); 

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
    // $ANTLR end "rule__ElkNodeEdges__Group_1__0__Impl"


    // $ANTLR start "rule__ElkNodeEdges__Group_1__1"
    // InternalElkGraphJson.g:4751:1: rule__ElkNodeEdges__Group_1__1 : rule__ElkNodeEdges__Group_1__1__Impl ;
    public final void rule__ElkNodeEdges__Group_1__1() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:4755:1: ( rule__ElkNodeEdges__Group_1__1__Impl )
            // InternalElkGraphJson.g:4756:2: rule__ElkNodeEdges__Group_1__1__Impl
            {
            pushFollow(FOLLOW_2);
            rule__ElkNodeEdges__Group_1__1__Impl();

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
    // $ANTLR end "rule__ElkNodeEdges__Group_1__1"


    // $ANTLR start "rule__ElkNodeEdges__Group_1__1__Impl"
    // InternalElkGraphJson.g:4762:1: rule__ElkNodeEdges__Group_1__1__Impl : ( ( rule__ElkNodeEdges__Group_1_1__0 )* ) ;
    public final void rule__ElkNodeEdges__Group_1__1__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:4766:1: ( ( ( rule__ElkNodeEdges__Group_1_1__0 )* ) )
            // InternalElkGraphJson.g:4767:1: ( ( rule__ElkNodeEdges__Group_1_1__0 )* )
            {
            // InternalElkGraphJson.g:4767:1: ( ( rule__ElkNodeEdges__Group_1_1__0 )* )
            // InternalElkGraphJson.g:4768:2: ( rule__ElkNodeEdges__Group_1_1__0 )*
            {
             before(grammarAccess.getElkNodeEdgesAccess().getGroup_1_1()); 
            // InternalElkGraphJson.g:4769:2: ( rule__ElkNodeEdges__Group_1_1__0 )*
            loop49:
            do {
                int alt49=2;
                int LA49_0 = input.LA(1);

                if ( (LA49_0==59) ) {
                    int LA49_1 = input.LA(2);

                    if ( (LA49_1==58) ) {
                        alt49=1;
                    }


                }


                switch (alt49) {
            	case 1 :
            	    // InternalElkGraphJson.g:4769:3: rule__ElkNodeEdges__Group_1_1__0
            	    {
            	    pushFollow(FOLLOW_6);
            	    rule__ElkNodeEdges__Group_1_1__0();

            	    state._fsp--;


            	    }
            	    break;

            	default :
            	    break loop49;
                }
            } while (true);

             after(grammarAccess.getElkNodeEdgesAccess().getGroup_1_1()); 

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
    // $ANTLR end "rule__ElkNodeEdges__Group_1__1__Impl"


    // $ANTLR start "rule__ElkNodeEdges__Group_1_1__0"
    // InternalElkGraphJson.g:4778:1: rule__ElkNodeEdges__Group_1_1__0 : rule__ElkNodeEdges__Group_1_1__0__Impl rule__ElkNodeEdges__Group_1_1__1 ;
    public final void rule__ElkNodeEdges__Group_1_1__0() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:4782:1: ( rule__ElkNodeEdges__Group_1_1__0__Impl rule__ElkNodeEdges__Group_1_1__1 )
            // InternalElkGraphJson.g:4783:2: rule__ElkNodeEdges__Group_1_1__0__Impl rule__ElkNodeEdges__Group_1_1__1
            {
            pushFollow(FOLLOW_3);
            rule__ElkNodeEdges__Group_1_1__0__Impl();

            state._fsp--;

            pushFollow(FOLLOW_2);
            rule__ElkNodeEdges__Group_1_1__1();

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
    // $ANTLR end "rule__ElkNodeEdges__Group_1_1__0"


    // $ANTLR start "rule__ElkNodeEdges__Group_1_1__0__Impl"
    // InternalElkGraphJson.g:4790:1: rule__ElkNodeEdges__Group_1_1__0__Impl : ( ',' ) ;
    public final void rule__ElkNodeEdges__Group_1_1__0__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:4794:1: ( ( ',' ) )
            // InternalElkGraphJson.g:4795:1: ( ',' )
            {
            // InternalElkGraphJson.g:4795:1: ( ',' )
            // InternalElkGraphJson.g:4796:2: ','
            {
             before(grammarAccess.getElkNodeEdgesAccess().getCommaKeyword_1_1_0()); 
            match(input,59,FOLLOW_2); 
             after(grammarAccess.getElkNodeEdgesAccess().getCommaKeyword_1_1_0()); 

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
    // $ANTLR end "rule__ElkNodeEdges__Group_1_1__0__Impl"


    // $ANTLR start "rule__ElkNodeEdges__Group_1_1__1"
    // InternalElkGraphJson.g:4805:1: rule__ElkNodeEdges__Group_1_1__1 : rule__ElkNodeEdges__Group_1_1__1__Impl ;
    public final void rule__ElkNodeEdges__Group_1_1__1() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:4809:1: ( rule__ElkNodeEdges__Group_1_1__1__Impl )
            // InternalElkGraphJson.g:4810:2: rule__ElkNodeEdges__Group_1_1__1__Impl
            {
            pushFollow(FOLLOW_2);
            rule__ElkNodeEdges__Group_1_1__1__Impl();

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
    // $ANTLR end "rule__ElkNodeEdges__Group_1_1__1"


    // $ANTLR start "rule__ElkNodeEdges__Group_1_1__1__Impl"
    // InternalElkGraphJson.g:4816:1: rule__ElkNodeEdges__Group_1_1__1__Impl : ( ( rule__ElkNodeEdges__ContainedEdgesAssignment_1_1_1 ) ) ;
    public final void rule__ElkNodeEdges__Group_1_1__1__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:4820:1: ( ( ( rule__ElkNodeEdges__ContainedEdgesAssignment_1_1_1 ) ) )
            // InternalElkGraphJson.g:4821:1: ( ( rule__ElkNodeEdges__ContainedEdgesAssignment_1_1_1 ) )
            {
            // InternalElkGraphJson.g:4821:1: ( ( rule__ElkNodeEdges__ContainedEdgesAssignment_1_1_1 ) )
            // InternalElkGraphJson.g:4822:2: ( rule__ElkNodeEdges__ContainedEdgesAssignment_1_1_1 )
            {
             before(grammarAccess.getElkNodeEdgesAccess().getContainedEdgesAssignment_1_1_1()); 
            // InternalElkGraphJson.g:4823:2: ( rule__ElkNodeEdges__ContainedEdgesAssignment_1_1_1 )
            // InternalElkGraphJson.g:4823:3: rule__ElkNodeEdges__ContainedEdgesAssignment_1_1_1
            {
            pushFollow(FOLLOW_2);
            rule__ElkNodeEdges__ContainedEdgesAssignment_1_1_1();

            state._fsp--;


            }

             after(grammarAccess.getElkNodeEdgesAccess().getContainedEdgesAssignment_1_1_1()); 

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
    // $ANTLR end "rule__ElkNodeEdges__Group_1_1__1__Impl"


    // $ANTLR start "rule__ElkGraphElementLabels__Group__0"
    // InternalElkGraphJson.g:4832:1: rule__ElkGraphElementLabels__Group__0 : rule__ElkGraphElementLabels__Group__0__Impl rule__ElkGraphElementLabels__Group__1 ;
    public final void rule__ElkGraphElementLabels__Group__0() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:4836:1: ( rule__ElkGraphElementLabels__Group__0__Impl rule__ElkGraphElementLabels__Group__1 )
            // InternalElkGraphJson.g:4837:2: rule__ElkGraphElementLabels__Group__0__Impl rule__ElkGraphElementLabels__Group__1
            {
            pushFollow(FOLLOW_15);
            rule__ElkGraphElementLabels__Group__0__Impl();

            state._fsp--;

            pushFollow(FOLLOW_2);
            rule__ElkGraphElementLabels__Group__1();

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
    // $ANTLR end "rule__ElkGraphElementLabels__Group__0"


    // $ANTLR start "rule__ElkGraphElementLabels__Group__0__Impl"
    // InternalElkGraphJson.g:4844:1: rule__ElkGraphElementLabels__Group__0__Impl : ( '[' ) ;
    public final void rule__ElkGraphElementLabels__Group__0__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:4848:1: ( ( '[' ) )
            // InternalElkGraphJson.g:4849:1: ( '[' )
            {
            // InternalElkGraphJson.g:4849:1: ( '[' )
            // InternalElkGraphJson.g:4850:2: '['
            {
             before(grammarAccess.getElkGraphElementLabelsAccess().getLeftSquareBracketKeyword_0()); 
            match(input,62,FOLLOW_2); 
             after(grammarAccess.getElkGraphElementLabelsAccess().getLeftSquareBracketKeyword_0()); 

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
    // $ANTLR end "rule__ElkGraphElementLabels__Group__0__Impl"


    // $ANTLR start "rule__ElkGraphElementLabels__Group__1"
    // InternalElkGraphJson.g:4859:1: rule__ElkGraphElementLabels__Group__1 : rule__ElkGraphElementLabels__Group__1__Impl rule__ElkGraphElementLabels__Group__2 ;
    public final void rule__ElkGraphElementLabels__Group__1() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:4863:1: ( rule__ElkGraphElementLabels__Group__1__Impl rule__ElkGraphElementLabels__Group__2 )
            // InternalElkGraphJson.g:4864:2: rule__ElkGraphElementLabels__Group__1__Impl rule__ElkGraphElementLabels__Group__2
            {
            pushFollow(FOLLOW_15);
            rule__ElkGraphElementLabels__Group__1__Impl();

            state._fsp--;

            pushFollow(FOLLOW_2);
            rule__ElkGraphElementLabels__Group__2();

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
    // $ANTLR end "rule__ElkGraphElementLabels__Group__1"


    // $ANTLR start "rule__ElkGraphElementLabels__Group__1__Impl"
    // InternalElkGraphJson.g:4871:1: rule__ElkGraphElementLabels__Group__1__Impl : ( ( rule__ElkGraphElementLabels__Group_1__0 )? ) ;
    public final void rule__ElkGraphElementLabels__Group__1__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:4875:1: ( ( ( rule__ElkGraphElementLabels__Group_1__0 )? ) )
            // InternalElkGraphJson.g:4876:1: ( ( rule__ElkGraphElementLabels__Group_1__0 )? )
            {
            // InternalElkGraphJson.g:4876:1: ( ( rule__ElkGraphElementLabels__Group_1__0 )? )
            // InternalElkGraphJson.g:4877:2: ( rule__ElkGraphElementLabels__Group_1__0 )?
            {
             before(grammarAccess.getElkGraphElementLabelsAccess().getGroup_1()); 
            // InternalElkGraphJson.g:4878:2: ( rule__ElkGraphElementLabels__Group_1__0 )?
            int alt50=2;
            int LA50_0 = input.LA(1);

            if ( (LA50_0==58) ) {
                alt50=1;
            }
            switch (alt50) {
                case 1 :
                    // InternalElkGraphJson.g:4878:3: rule__ElkGraphElementLabels__Group_1__0
                    {
                    pushFollow(FOLLOW_2);
                    rule__ElkGraphElementLabels__Group_1__0();

                    state._fsp--;


                    }
                    break;

            }

             after(grammarAccess.getElkGraphElementLabelsAccess().getGroup_1()); 

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
    // $ANTLR end "rule__ElkGraphElementLabels__Group__1__Impl"


    // $ANTLR start "rule__ElkGraphElementLabels__Group__2"
    // InternalElkGraphJson.g:4886:1: rule__ElkGraphElementLabels__Group__2 : rule__ElkGraphElementLabels__Group__2__Impl rule__ElkGraphElementLabels__Group__3 ;
    public final void rule__ElkGraphElementLabels__Group__2() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:4890:1: ( rule__ElkGraphElementLabels__Group__2__Impl rule__ElkGraphElementLabels__Group__3 )
            // InternalElkGraphJson.g:4891:2: rule__ElkGraphElementLabels__Group__2__Impl rule__ElkGraphElementLabels__Group__3
            {
            pushFollow(FOLLOW_15);
            rule__ElkGraphElementLabels__Group__2__Impl();

            state._fsp--;

            pushFollow(FOLLOW_2);
            rule__ElkGraphElementLabels__Group__3();

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
    // $ANTLR end "rule__ElkGraphElementLabels__Group__2"


    // $ANTLR start "rule__ElkGraphElementLabels__Group__2__Impl"
    // InternalElkGraphJson.g:4898:1: rule__ElkGraphElementLabels__Group__2__Impl : ( ( ',' )? ) ;
    public final void rule__ElkGraphElementLabels__Group__2__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:4902:1: ( ( ( ',' )? ) )
            // InternalElkGraphJson.g:4903:1: ( ( ',' )? )
            {
            // InternalElkGraphJson.g:4903:1: ( ( ',' )? )
            // InternalElkGraphJson.g:4904:2: ( ',' )?
            {
             before(grammarAccess.getElkGraphElementLabelsAccess().getCommaKeyword_2()); 
            // InternalElkGraphJson.g:4905:2: ( ',' )?
            int alt51=2;
            int LA51_0 = input.LA(1);

            if ( (LA51_0==59) ) {
                alt51=1;
            }
            switch (alt51) {
                case 1 :
                    // InternalElkGraphJson.g:4905:3: ','
                    {
                    match(input,59,FOLLOW_2); 

                    }
                    break;

            }

             after(grammarAccess.getElkGraphElementLabelsAccess().getCommaKeyword_2()); 

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
    // $ANTLR end "rule__ElkGraphElementLabels__Group__2__Impl"


    // $ANTLR start "rule__ElkGraphElementLabels__Group__3"
    // InternalElkGraphJson.g:4913:1: rule__ElkGraphElementLabels__Group__3 : rule__ElkGraphElementLabels__Group__3__Impl ;
    public final void rule__ElkGraphElementLabels__Group__3() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:4917:1: ( rule__ElkGraphElementLabels__Group__3__Impl )
            // InternalElkGraphJson.g:4918:2: rule__ElkGraphElementLabels__Group__3__Impl
            {
            pushFollow(FOLLOW_2);
            rule__ElkGraphElementLabels__Group__3__Impl();

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
    // $ANTLR end "rule__ElkGraphElementLabels__Group__3"


    // $ANTLR start "rule__ElkGraphElementLabels__Group__3__Impl"
    // InternalElkGraphJson.g:4924:1: rule__ElkGraphElementLabels__Group__3__Impl : ( ']' ) ;
    public final void rule__ElkGraphElementLabels__Group__3__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:4928:1: ( ( ']' ) )
            // InternalElkGraphJson.g:4929:1: ( ']' )
            {
            // InternalElkGraphJson.g:4929:1: ( ']' )
            // InternalElkGraphJson.g:4930:2: ']'
            {
             before(grammarAccess.getElkGraphElementLabelsAccess().getRightSquareBracketKeyword_3()); 
            match(input,63,FOLLOW_2); 
             after(grammarAccess.getElkGraphElementLabelsAccess().getRightSquareBracketKeyword_3()); 

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
    // $ANTLR end "rule__ElkGraphElementLabels__Group__3__Impl"


    // $ANTLR start "rule__ElkGraphElementLabels__Group_1__0"
    // InternalElkGraphJson.g:4940:1: rule__ElkGraphElementLabels__Group_1__0 : rule__ElkGraphElementLabels__Group_1__0__Impl rule__ElkGraphElementLabels__Group_1__1 ;
    public final void rule__ElkGraphElementLabels__Group_1__0() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:4944:1: ( rule__ElkGraphElementLabels__Group_1__0__Impl rule__ElkGraphElementLabels__Group_1__1 )
            // InternalElkGraphJson.g:4945:2: rule__ElkGraphElementLabels__Group_1__0__Impl rule__ElkGraphElementLabels__Group_1__1
            {
            pushFollow(FOLLOW_5);
            rule__ElkGraphElementLabels__Group_1__0__Impl();

            state._fsp--;

            pushFollow(FOLLOW_2);
            rule__ElkGraphElementLabels__Group_1__1();

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
    // $ANTLR end "rule__ElkGraphElementLabels__Group_1__0"


    // $ANTLR start "rule__ElkGraphElementLabels__Group_1__0__Impl"
    // InternalElkGraphJson.g:4952:1: rule__ElkGraphElementLabels__Group_1__0__Impl : ( ( rule__ElkGraphElementLabels__LabelsAssignment_1_0 ) ) ;
    public final void rule__ElkGraphElementLabels__Group_1__0__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:4956:1: ( ( ( rule__ElkGraphElementLabels__LabelsAssignment_1_0 ) ) )
            // InternalElkGraphJson.g:4957:1: ( ( rule__ElkGraphElementLabels__LabelsAssignment_1_0 ) )
            {
            // InternalElkGraphJson.g:4957:1: ( ( rule__ElkGraphElementLabels__LabelsAssignment_1_0 ) )
            // InternalElkGraphJson.g:4958:2: ( rule__ElkGraphElementLabels__LabelsAssignment_1_0 )
            {
             before(grammarAccess.getElkGraphElementLabelsAccess().getLabelsAssignment_1_0()); 
            // InternalElkGraphJson.g:4959:2: ( rule__ElkGraphElementLabels__LabelsAssignment_1_0 )
            // InternalElkGraphJson.g:4959:3: rule__ElkGraphElementLabels__LabelsAssignment_1_0
            {
            pushFollow(FOLLOW_2);
            rule__ElkGraphElementLabels__LabelsAssignment_1_0();

            state._fsp--;


            }

             after(grammarAccess.getElkGraphElementLabelsAccess().getLabelsAssignment_1_0()); 

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
    // $ANTLR end "rule__ElkGraphElementLabels__Group_1__0__Impl"


    // $ANTLR start "rule__ElkGraphElementLabels__Group_1__1"
    // InternalElkGraphJson.g:4967:1: rule__ElkGraphElementLabels__Group_1__1 : rule__ElkGraphElementLabels__Group_1__1__Impl ;
    public final void rule__ElkGraphElementLabels__Group_1__1() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:4971:1: ( rule__ElkGraphElementLabels__Group_1__1__Impl )
            // InternalElkGraphJson.g:4972:2: rule__ElkGraphElementLabels__Group_1__1__Impl
            {
            pushFollow(FOLLOW_2);
            rule__ElkGraphElementLabels__Group_1__1__Impl();

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
    // $ANTLR end "rule__ElkGraphElementLabels__Group_1__1"


    // $ANTLR start "rule__ElkGraphElementLabels__Group_1__1__Impl"
    // InternalElkGraphJson.g:4978:1: rule__ElkGraphElementLabels__Group_1__1__Impl : ( ( rule__ElkGraphElementLabels__Group_1_1__0 )* ) ;
    public final void rule__ElkGraphElementLabels__Group_1__1__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:4982:1: ( ( ( rule__ElkGraphElementLabels__Group_1_1__0 )* ) )
            // InternalElkGraphJson.g:4983:1: ( ( rule__ElkGraphElementLabels__Group_1_1__0 )* )
            {
            // InternalElkGraphJson.g:4983:1: ( ( rule__ElkGraphElementLabels__Group_1_1__0 )* )
            // InternalElkGraphJson.g:4984:2: ( rule__ElkGraphElementLabels__Group_1_1__0 )*
            {
             before(grammarAccess.getElkGraphElementLabelsAccess().getGroup_1_1()); 
            // InternalElkGraphJson.g:4985:2: ( rule__ElkGraphElementLabels__Group_1_1__0 )*
            loop52:
            do {
                int alt52=2;
                int LA52_0 = input.LA(1);

                if ( (LA52_0==59) ) {
                    int LA52_1 = input.LA(2);

                    if ( (LA52_1==58) ) {
                        alt52=1;
                    }


                }


                switch (alt52) {
            	case 1 :
            	    // InternalElkGraphJson.g:4985:3: rule__ElkGraphElementLabels__Group_1_1__0
            	    {
            	    pushFollow(FOLLOW_6);
            	    rule__ElkGraphElementLabels__Group_1_1__0();

            	    state._fsp--;


            	    }
            	    break;

            	default :
            	    break loop52;
                }
            } while (true);

             after(grammarAccess.getElkGraphElementLabelsAccess().getGroup_1_1()); 

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
    // $ANTLR end "rule__ElkGraphElementLabels__Group_1__1__Impl"


    // $ANTLR start "rule__ElkGraphElementLabels__Group_1_1__0"
    // InternalElkGraphJson.g:4994:1: rule__ElkGraphElementLabels__Group_1_1__0 : rule__ElkGraphElementLabels__Group_1_1__0__Impl rule__ElkGraphElementLabels__Group_1_1__1 ;
    public final void rule__ElkGraphElementLabels__Group_1_1__0() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:4998:1: ( rule__ElkGraphElementLabels__Group_1_1__0__Impl rule__ElkGraphElementLabels__Group_1_1__1 )
            // InternalElkGraphJson.g:4999:2: rule__ElkGraphElementLabels__Group_1_1__0__Impl rule__ElkGraphElementLabels__Group_1_1__1
            {
            pushFollow(FOLLOW_3);
            rule__ElkGraphElementLabels__Group_1_1__0__Impl();

            state._fsp--;

            pushFollow(FOLLOW_2);
            rule__ElkGraphElementLabels__Group_1_1__1();

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
    // $ANTLR end "rule__ElkGraphElementLabels__Group_1_1__0"


    // $ANTLR start "rule__ElkGraphElementLabels__Group_1_1__0__Impl"
    // InternalElkGraphJson.g:5006:1: rule__ElkGraphElementLabels__Group_1_1__0__Impl : ( ',' ) ;
    public final void rule__ElkGraphElementLabels__Group_1_1__0__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:5010:1: ( ( ',' ) )
            // InternalElkGraphJson.g:5011:1: ( ',' )
            {
            // InternalElkGraphJson.g:5011:1: ( ',' )
            // InternalElkGraphJson.g:5012:2: ','
            {
             before(grammarAccess.getElkGraphElementLabelsAccess().getCommaKeyword_1_1_0()); 
            match(input,59,FOLLOW_2); 
             after(grammarAccess.getElkGraphElementLabelsAccess().getCommaKeyword_1_1_0()); 

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
    // $ANTLR end "rule__ElkGraphElementLabels__Group_1_1__0__Impl"


    // $ANTLR start "rule__ElkGraphElementLabels__Group_1_1__1"
    // InternalElkGraphJson.g:5021:1: rule__ElkGraphElementLabels__Group_1_1__1 : rule__ElkGraphElementLabels__Group_1_1__1__Impl ;
    public final void rule__ElkGraphElementLabels__Group_1_1__1() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:5025:1: ( rule__ElkGraphElementLabels__Group_1_1__1__Impl )
            // InternalElkGraphJson.g:5026:2: rule__ElkGraphElementLabels__Group_1_1__1__Impl
            {
            pushFollow(FOLLOW_2);
            rule__ElkGraphElementLabels__Group_1_1__1__Impl();

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
    // $ANTLR end "rule__ElkGraphElementLabels__Group_1_1__1"


    // $ANTLR start "rule__ElkGraphElementLabels__Group_1_1__1__Impl"
    // InternalElkGraphJson.g:5032:1: rule__ElkGraphElementLabels__Group_1_1__1__Impl : ( ( rule__ElkGraphElementLabels__LabelsAssignment_1_1_1 ) ) ;
    public final void rule__ElkGraphElementLabels__Group_1_1__1__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:5036:1: ( ( ( rule__ElkGraphElementLabels__LabelsAssignment_1_1_1 ) ) )
            // InternalElkGraphJson.g:5037:1: ( ( rule__ElkGraphElementLabels__LabelsAssignment_1_1_1 ) )
            {
            // InternalElkGraphJson.g:5037:1: ( ( rule__ElkGraphElementLabels__LabelsAssignment_1_1_1 ) )
            // InternalElkGraphJson.g:5038:2: ( rule__ElkGraphElementLabels__LabelsAssignment_1_1_1 )
            {
             before(grammarAccess.getElkGraphElementLabelsAccess().getLabelsAssignment_1_1_1()); 
            // InternalElkGraphJson.g:5039:2: ( rule__ElkGraphElementLabels__LabelsAssignment_1_1_1 )
            // InternalElkGraphJson.g:5039:3: rule__ElkGraphElementLabels__LabelsAssignment_1_1_1
            {
            pushFollow(FOLLOW_2);
            rule__ElkGraphElementLabels__LabelsAssignment_1_1_1();

            state._fsp--;


            }

             after(grammarAccess.getElkGraphElementLabelsAccess().getLabelsAssignment_1_1_1()); 

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
    // $ANTLR end "rule__ElkGraphElementLabels__Group_1_1__1__Impl"


    // $ANTLR start "rule__ElkGraphElementProperties__Group__0"
    // InternalElkGraphJson.g:5048:1: rule__ElkGraphElementProperties__Group__0 : rule__ElkGraphElementProperties__Group__0__Impl rule__ElkGraphElementProperties__Group__1 ;
    public final void rule__ElkGraphElementProperties__Group__0() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:5052:1: ( rule__ElkGraphElementProperties__Group__0__Impl rule__ElkGraphElementProperties__Group__1 )
            // InternalElkGraphJson.g:5053:2: rule__ElkGraphElementProperties__Group__0__Impl rule__ElkGraphElementProperties__Group__1
            {
            pushFollow(FOLLOW_16);
            rule__ElkGraphElementProperties__Group__0__Impl();

            state._fsp--;

            pushFollow(FOLLOW_2);
            rule__ElkGraphElementProperties__Group__1();

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
    // $ANTLR end "rule__ElkGraphElementProperties__Group__0"


    // $ANTLR start "rule__ElkGraphElementProperties__Group__0__Impl"
    // InternalElkGraphJson.g:5060:1: rule__ElkGraphElementProperties__Group__0__Impl : ( '{' ) ;
    public final void rule__ElkGraphElementProperties__Group__0__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:5064:1: ( ( '{' ) )
            // InternalElkGraphJson.g:5065:1: ( '{' )
            {
            // InternalElkGraphJson.g:5065:1: ( '{' )
            // InternalElkGraphJson.g:5066:2: '{'
            {
             before(grammarAccess.getElkGraphElementPropertiesAccess().getLeftCurlyBracketKeyword_0()); 
            match(input,58,FOLLOW_2); 
             after(grammarAccess.getElkGraphElementPropertiesAccess().getLeftCurlyBracketKeyword_0()); 

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
    // $ANTLR end "rule__ElkGraphElementProperties__Group__0__Impl"


    // $ANTLR start "rule__ElkGraphElementProperties__Group__1"
    // InternalElkGraphJson.g:5075:1: rule__ElkGraphElementProperties__Group__1 : rule__ElkGraphElementProperties__Group__1__Impl rule__ElkGraphElementProperties__Group__2 ;
    public final void rule__ElkGraphElementProperties__Group__1() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:5079:1: ( rule__ElkGraphElementProperties__Group__1__Impl rule__ElkGraphElementProperties__Group__2 )
            // InternalElkGraphJson.g:5080:2: rule__ElkGraphElementProperties__Group__1__Impl rule__ElkGraphElementProperties__Group__2
            {
            pushFollow(FOLLOW_16);
            rule__ElkGraphElementProperties__Group__1__Impl();

            state._fsp--;

            pushFollow(FOLLOW_2);
            rule__ElkGraphElementProperties__Group__2();

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
    // $ANTLR end "rule__ElkGraphElementProperties__Group__1"


    // $ANTLR start "rule__ElkGraphElementProperties__Group__1__Impl"
    // InternalElkGraphJson.g:5087:1: rule__ElkGraphElementProperties__Group__1__Impl : ( ( rule__ElkGraphElementProperties__Group_1__0 )? ) ;
    public final void rule__ElkGraphElementProperties__Group__1__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:5091:1: ( ( ( rule__ElkGraphElementProperties__Group_1__0 )? ) )
            // InternalElkGraphJson.g:5092:1: ( ( rule__ElkGraphElementProperties__Group_1__0 )? )
            {
            // InternalElkGraphJson.g:5092:1: ( ( rule__ElkGraphElementProperties__Group_1__0 )? )
            // InternalElkGraphJson.g:5093:2: ( rule__ElkGraphElementProperties__Group_1__0 )?
            {
             before(grammarAccess.getElkGraphElementPropertiesAccess().getGroup_1()); 
            // InternalElkGraphJson.g:5094:2: ( rule__ElkGraphElementProperties__Group_1__0 )?
            int alt53=2;
            int LA53_0 = input.LA(1);

            if ( ((LA53_0>=RULE_STRING && LA53_0<=RULE_ID)) ) {
                alt53=1;
            }
            switch (alt53) {
                case 1 :
                    // InternalElkGraphJson.g:5094:3: rule__ElkGraphElementProperties__Group_1__0
                    {
                    pushFollow(FOLLOW_2);
                    rule__ElkGraphElementProperties__Group_1__0();

                    state._fsp--;


                    }
                    break;

            }

             after(grammarAccess.getElkGraphElementPropertiesAccess().getGroup_1()); 

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
    // $ANTLR end "rule__ElkGraphElementProperties__Group__1__Impl"


    // $ANTLR start "rule__ElkGraphElementProperties__Group__2"
    // InternalElkGraphJson.g:5102:1: rule__ElkGraphElementProperties__Group__2 : rule__ElkGraphElementProperties__Group__2__Impl rule__ElkGraphElementProperties__Group__3 ;
    public final void rule__ElkGraphElementProperties__Group__2() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:5106:1: ( rule__ElkGraphElementProperties__Group__2__Impl rule__ElkGraphElementProperties__Group__3 )
            // InternalElkGraphJson.g:5107:2: rule__ElkGraphElementProperties__Group__2__Impl rule__ElkGraphElementProperties__Group__3
            {
            pushFollow(FOLLOW_16);
            rule__ElkGraphElementProperties__Group__2__Impl();

            state._fsp--;

            pushFollow(FOLLOW_2);
            rule__ElkGraphElementProperties__Group__3();

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
    // $ANTLR end "rule__ElkGraphElementProperties__Group__2"


    // $ANTLR start "rule__ElkGraphElementProperties__Group__2__Impl"
    // InternalElkGraphJson.g:5114:1: rule__ElkGraphElementProperties__Group__2__Impl : ( ( ',' )? ) ;
    public final void rule__ElkGraphElementProperties__Group__2__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:5118:1: ( ( ( ',' )? ) )
            // InternalElkGraphJson.g:5119:1: ( ( ',' )? )
            {
            // InternalElkGraphJson.g:5119:1: ( ( ',' )? )
            // InternalElkGraphJson.g:5120:2: ( ',' )?
            {
             before(grammarAccess.getElkGraphElementPropertiesAccess().getCommaKeyword_2()); 
            // InternalElkGraphJson.g:5121:2: ( ',' )?
            int alt54=2;
            int LA54_0 = input.LA(1);

            if ( (LA54_0==59) ) {
                alt54=1;
            }
            switch (alt54) {
                case 1 :
                    // InternalElkGraphJson.g:5121:3: ','
                    {
                    match(input,59,FOLLOW_2); 

                    }
                    break;

            }

             after(grammarAccess.getElkGraphElementPropertiesAccess().getCommaKeyword_2()); 

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
    // $ANTLR end "rule__ElkGraphElementProperties__Group__2__Impl"


    // $ANTLR start "rule__ElkGraphElementProperties__Group__3"
    // InternalElkGraphJson.g:5129:1: rule__ElkGraphElementProperties__Group__3 : rule__ElkGraphElementProperties__Group__3__Impl ;
    public final void rule__ElkGraphElementProperties__Group__3() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:5133:1: ( rule__ElkGraphElementProperties__Group__3__Impl )
            // InternalElkGraphJson.g:5134:2: rule__ElkGraphElementProperties__Group__3__Impl
            {
            pushFollow(FOLLOW_2);
            rule__ElkGraphElementProperties__Group__3__Impl();

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
    // $ANTLR end "rule__ElkGraphElementProperties__Group__3"


    // $ANTLR start "rule__ElkGraphElementProperties__Group__3__Impl"
    // InternalElkGraphJson.g:5140:1: rule__ElkGraphElementProperties__Group__3__Impl : ( '}' ) ;
    public final void rule__ElkGraphElementProperties__Group__3__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:5144:1: ( ( '}' ) )
            // InternalElkGraphJson.g:5145:1: ( '}' )
            {
            // InternalElkGraphJson.g:5145:1: ( '}' )
            // InternalElkGraphJson.g:5146:2: '}'
            {
             before(grammarAccess.getElkGraphElementPropertiesAccess().getRightCurlyBracketKeyword_3()); 
            match(input,60,FOLLOW_2); 
             after(grammarAccess.getElkGraphElementPropertiesAccess().getRightCurlyBracketKeyword_3()); 

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
    // $ANTLR end "rule__ElkGraphElementProperties__Group__3__Impl"


    // $ANTLR start "rule__ElkGraphElementProperties__Group_1__0"
    // InternalElkGraphJson.g:5156:1: rule__ElkGraphElementProperties__Group_1__0 : rule__ElkGraphElementProperties__Group_1__0__Impl rule__ElkGraphElementProperties__Group_1__1 ;
    public final void rule__ElkGraphElementProperties__Group_1__0() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:5160:1: ( rule__ElkGraphElementProperties__Group_1__0__Impl rule__ElkGraphElementProperties__Group_1__1 )
            // InternalElkGraphJson.g:5161:2: rule__ElkGraphElementProperties__Group_1__0__Impl rule__ElkGraphElementProperties__Group_1__1
            {
            pushFollow(FOLLOW_5);
            rule__ElkGraphElementProperties__Group_1__0__Impl();

            state._fsp--;

            pushFollow(FOLLOW_2);
            rule__ElkGraphElementProperties__Group_1__1();

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
    // $ANTLR end "rule__ElkGraphElementProperties__Group_1__0"


    // $ANTLR start "rule__ElkGraphElementProperties__Group_1__0__Impl"
    // InternalElkGraphJson.g:5168:1: rule__ElkGraphElementProperties__Group_1__0__Impl : ( ( rule__ElkGraphElementProperties__PropertiesAssignment_1_0 ) ) ;
    public final void rule__ElkGraphElementProperties__Group_1__0__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:5172:1: ( ( ( rule__ElkGraphElementProperties__PropertiesAssignment_1_0 ) ) )
            // InternalElkGraphJson.g:5173:1: ( ( rule__ElkGraphElementProperties__PropertiesAssignment_1_0 ) )
            {
            // InternalElkGraphJson.g:5173:1: ( ( rule__ElkGraphElementProperties__PropertiesAssignment_1_0 ) )
            // InternalElkGraphJson.g:5174:2: ( rule__ElkGraphElementProperties__PropertiesAssignment_1_0 )
            {
             before(grammarAccess.getElkGraphElementPropertiesAccess().getPropertiesAssignment_1_0()); 
            // InternalElkGraphJson.g:5175:2: ( rule__ElkGraphElementProperties__PropertiesAssignment_1_0 )
            // InternalElkGraphJson.g:5175:3: rule__ElkGraphElementProperties__PropertiesAssignment_1_0
            {
            pushFollow(FOLLOW_2);
            rule__ElkGraphElementProperties__PropertiesAssignment_1_0();

            state._fsp--;


            }

             after(grammarAccess.getElkGraphElementPropertiesAccess().getPropertiesAssignment_1_0()); 

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
    // $ANTLR end "rule__ElkGraphElementProperties__Group_1__0__Impl"


    // $ANTLR start "rule__ElkGraphElementProperties__Group_1__1"
    // InternalElkGraphJson.g:5183:1: rule__ElkGraphElementProperties__Group_1__1 : rule__ElkGraphElementProperties__Group_1__1__Impl ;
    public final void rule__ElkGraphElementProperties__Group_1__1() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:5187:1: ( rule__ElkGraphElementProperties__Group_1__1__Impl )
            // InternalElkGraphJson.g:5188:2: rule__ElkGraphElementProperties__Group_1__1__Impl
            {
            pushFollow(FOLLOW_2);
            rule__ElkGraphElementProperties__Group_1__1__Impl();

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
    // $ANTLR end "rule__ElkGraphElementProperties__Group_1__1"


    // $ANTLR start "rule__ElkGraphElementProperties__Group_1__1__Impl"
    // InternalElkGraphJson.g:5194:1: rule__ElkGraphElementProperties__Group_1__1__Impl : ( ( rule__ElkGraphElementProperties__Group_1_1__0 )* ) ;
    public final void rule__ElkGraphElementProperties__Group_1__1__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:5198:1: ( ( ( rule__ElkGraphElementProperties__Group_1_1__0 )* ) )
            // InternalElkGraphJson.g:5199:1: ( ( rule__ElkGraphElementProperties__Group_1_1__0 )* )
            {
            // InternalElkGraphJson.g:5199:1: ( ( rule__ElkGraphElementProperties__Group_1_1__0 )* )
            // InternalElkGraphJson.g:5200:2: ( rule__ElkGraphElementProperties__Group_1_1__0 )*
            {
             before(grammarAccess.getElkGraphElementPropertiesAccess().getGroup_1_1()); 
            // InternalElkGraphJson.g:5201:2: ( rule__ElkGraphElementProperties__Group_1_1__0 )*
            loop55:
            do {
                int alt55=2;
                int LA55_0 = input.LA(1);

                if ( (LA55_0==59) ) {
                    int LA55_1 = input.LA(2);

                    if ( ((LA55_1>=RULE_STRING && LA55_1<=RULE_ID)) ) {
                        alt55=1;
                    }


                }


                switch (alt55) {
            	case 1 :
            	    // InternalElkGraphJson.g:5201:3: rule__ElkGraphElementProperties__Group_1_1__0
            	    {
            	    pushFollow(FOLLOW_6);
            	    rule__ElkGraphElementProperties__Group_1_1__0();

            	    state._fsp--;


            	    }
            	    break;

            	default :
            	    break loop55;
                }
            } while (true);

             after(grammarAccess.getElkGraphElementPropertiesAccess().getGroup_1_1()); 

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
    // $ANTLR end "rule__ElkGraphElementProperties__Group_1__1__Impl"


    // $ANTLR start "rule__ElkGraphElementProperties__Group_1_1__0"
    // InternalElkGraphJson.g:5210:1: rule__ElkGraphElementProperties__Group_1_1__0 : rule__ElkGraphElementProperties__Group_1_1__0__Impl rule__ElkGraphElementProperties__Group_1_1__1 ;
    public final void rule__ElkGraphElementProperties__Group_1_1__0() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:5214:1: ( rule__ElkGraphElementProperties__Group_1_1__0__Impl rule__ElkGraphElementProperties__Group_1_1__1 )
            // InternalElkGraphJson.g:5215:2: rule__ElkGraphElementProperties__Group_1_1__0__Impl rule__ElkGraphElementProperties__Group_1_1__1
            {
            pushFollow(FOLLOW_17);
            rule__ElkGraphElementProperties__Group_1_1__0__Impl();

            state._fsp--;

            pushFollow(FOLLOW_2);
            rule__ElkGraphElementProperties__Group_1_1__1();

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
    // $ANTLR end "rule__ElkGraphElementProperties__Group_1_1__0"


    // $ANTLR start "rule__ElkGraphElementProperties__Group_1_1__0__Impl"
    // InternalElkGraphJson.g:5222:1: rule__ElkGraphElementProperties__Group_1_1__0__Impl : ( ',' ) ;
    public final void rule__ElkGraphElementProperties__Group_1_1__0__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:5226:1: ( ( ',' ) )
            // InternalElkGraphJson.g:5227:1: ( ',' )
            {
            // InternalElkGraphJson.g:5227:1: ( ',' )
            // InternalElkGraphJson.g:5228:2: ','
            {
             before(grammarAccess.getElkGraphElementPropertiesAccess().getCommaKeyword_1_1_0()); 
            match(input,59,FOLLOW_2); 
             after(grammarAccess.getElkGraphElementPropertiesAccess().getCommaKeyword_1_1_0()); 

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
    // $ANTLR end "rule__ElkGraphElementProperties__Group_1_1__0__Impl"


    // $ANTLR start "rule__ElkGraphElementProperties__Group_1_1__1"
    // InternalElkGraphJson.g:5237:1: rule__ElkGraphElementProperties__Group_1_1__1 : rule__ElkGraphElementProperties__Group_1_1__1__Impl ;
    public final void rule__ElkGraphElementProperties__Group_1_1__1() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:5241:1: ( rule__ElkGraphElementProperties__Group_1_1__1__Impl )
            // InternalElkGraphJson.g:5242:2: rule__ElkGraphElementProperties__Group_1_1__1__Impl
            {
            pushFollow(FOLLOW_2);
            rule__ElkGraphElementProperties__Group_1_1__1__Impl();

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
    // $ANTLR end "rule__ElkGraphElementProperties__Group_1_1__1"


    // $ANTLR start "rule__ElkGraphElementProperties__Group_1_1__1__Impl"
    // InternalElkGraphJson.g:5248:1: rule__ElkGraphElementProperties__Group_1_1__1__Impl : ( ( rule__ElkGraphElementProperties__PropertiesAssignment_1_1_1 ) ) ;
    public final void rule__ElkGraphElementProperties__Group_1_1__1__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:5252:1: ( ( ( rule__ElkGraphElementProperties__PropertiesAssignment_1_1_1 ) ) )
            // InternalElkGraphJson.g:5253:1: ( ( rule__ElkGraphElementProperties__PropertiesAssignment_1_1_1 ) )
            {
            // InternalElkGraphJson.g:5253:1: ( ( rule__ElkGraphElementProperties__PropertiesAssignment_1_1_1 ) )
            // InternalElkGraphJson.g:5254:2: ( rule__ElkGraphElementProperties__PropertiesAssignment_1_1_1 )
            {
             before(grammarAccess.getElkGraphElementPropertiesAccess().getPropertiesAssignment_1_1_1()); 
            // InternalElkGraphJson.g:5255:2: ( rule__ElkGraphElementProperties__PropertiesAssignment_1_1_1 )
            // InternalElkGraphJson.g:5255:3: rule__ElkGraphElementProperties__PropertiesAssignment_1_1_1
            {
            pushFollow(FOLLOW_2);
            rule__ElkGraphElementProperties__PropertiesAssignment_1_1_1();

            state._fsp--;


            }

             after(grammarAccess.getElkGraphElementPropertiesAccess().getPropertiesAssignment_1_1_1()); 

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
    // $ANTLR end "rule__ElkGraphElementProperties__Group_1_1__1__Impl"


    // $ANTLR start "rule__ShapeElement__Group_0__0"
    // InternalElkGraphJson.g:5264:1: rule__ShapeElement__Group_0__0 : rule__ShapeElement__Group_0__0__Impl rule__ShapeElement__Group_0__1 ;
    public final void rule__ShapeElement__Group_0__0() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:5268:1: ( rule__ShapeElement__Group_0__0__Impl rule__ShapeElement__Group_0__1 )
            // InternalElkGraphJson.g:5269:2: rule__ShapeElement__Group_0__0__Impl rule__ShapeElement__Group_0__1
            {
            pushFollow(FOLLOW_8);
            rule__ShapeElement__Group_0__0__Impl();

            state._fsp--;

            pushFollow(FOLLOW_2);
            rule__ShapeElement__Group_0__1();

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
    // $ANTLR end "rule__ShapeElement__Group_0__0"


    // $ANTLR start "rule__ShapeElement__Group_0__0__Impl"
    // InternalElkGraphJson.g:5276:1: rule__ShapeElement__Group_0__0__Impl : ( ruleKeyX ) ;
    public final void rule__ShapeElement__Group_0__0__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:5280:1: ( ( ruleKeyX ) )
            // InternalElkGraphJson.g:5281:1: ( ruleKeyX )
            {
            // InternalElkGraphJson.g:5281:1: ( ruleKeyX )
            // InternalElkGraphJson.g:5282:2: ruleKeyX
            {
             before(grammarAccess.getShapeElementAccess().getKeyXParserRuleCall_0_0()); 
            pushFollow(FOLLOW_2);
            ruleKeyX();

            state._fsp--;

             after(grammarAccess.getShapeElementAccess().getKeyXParserRuleCall_0_0()); 

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
    // $ANTLR end "rule__ShapeElement__Group_0__0__Impl"


    // $ANTLR start "rule__ShapeElement__Group_0__1"
    // InternalElkGraphJson.g:5291:1: rule__ShapeElement__Group_0__1 : rule__ShapeElement__Group_0__1__Impl rule__ShapeElement__Group_0__2 ;
    public final void rule__ShapeElement__Group_0__1() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:5295:1: ( rule__ShapeElement__Group_0__1__Impl rule__ShapeElement__Group_0__2 )
            // InternalElkGraphJson.g:5296:2: rule__ShapeElement__Group_0__1__Impl rule__ShapeElement__Group_0__2
            {
            pushFollow(FOLLOW_18);
            rule__ShapeElement__Group_0__1__Impl();

            state._fsp--;

            pushFollow(FOLLOW_2);
            rule__ShapeElement__Group_0__2();

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
    // $ANTLR end "rule__ShapeElement__Group_0__1"


    // $ANTLR start "rule__ShapeElement__Group_0__1__Impl"
    // InternalElkGraphJson.g:5303:1: rule__ShapeElement__Group_0__1__Impl : ( ':' ) ;
    public final void rule__ShapeElement__Group_0__1__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:5307:1: ( ( ':' ) )
            // InternalElkGraphJson.g:5308:1: ( ':' )
            {
            // InternalElkGraphJson.g:5308:1: ( ':' )
            // InternalElkGraphJson.g:5309:2: ':'
            {
             before(grammarAccess.getShapeElementAccess().getColonKeyword_0_1()); 
            match(input,61,FOLLOW_2); 
             after(grammarAccess.getShapeElementAccess().getColonKeyword_0_1()); 

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
    // $ANTLR end "rule__ShapeElement__Group_0__1__Impl"


    // $ANTLR start "rule__ShapeElement__Group_0__2"
    // InternalElkGraphJson.g:5318:1: rule__ShapeElement__Group_0__2 : rule__ShapeElement__Group_0__2__Impl ;
    public final void rule__ShapeElement__Group_0__2() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:5322:1: ( rule__ShapeElement__Group_0__2__Impl )
            // InternalElkGraphJson.g:5323:2: rule__ShapeElement__Group_0__2__Impl
            {
            pushFollow(FOLLOW_2);
            rule__ShapeElement__Group_0__2__Impl();

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
    // $ANTLR end "rule__ShapeElement__Group_0__2"


    // $ANTLR start "rule__ShapeElement__Group_0__2__Impl"
    // InternalElkGraphJson.g:5329:1: rule__ShapeElement__Group_0__2__Impl : ( ( rule__ShapeElement__XAssignment_0_2 ) ) ;
    public final void rule__ShapeElement__Group_0__2__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:5333:1: ( ( ( rule__ShapeElement__XAssignment_0_2 ) ) )
            // InternalElkGraphJson.g:5334:1: ( ( rule__ShapeElement__XAssignment_0_2 ) )
            {
            // InternalElkGraphJson.g:5334:1: ( ( rule__ShapeElement__XAssignment_0_2 ) )
            // InternalElkGraphJson.g:5335:2: ( rule__ShapeElement__XAssignment_0_2 )
            {
             before(grammarAccess.getShapeElementAccess().getXAssignment_0_2()); 
            // InternalElkGraphJson.g:5336:2: ( rule__ShapeElement__XAssignment_0_2 )
            // InternalElkGraphJson.g:5336:3: rule__ShapeElement__XAssignment_0_2
            {
            pushFollow(FOLLOW_2);
            rule__ShapeElement__XAssignment_0_2();

            state._fsp--;


            }

             after(grammarAccess.getShapeElementAccess().getXAssignment_0_2()); 

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
    // $ANTLR end "rule__ShapeElement__Group_0__2__Impl"


    // $ANTLR start "rule__ShapeElement__Group_1__0"
    // InternalElkGraphJson.g:5345:1: rule__ShapeElement__Group_1__0 : rule__ShapeElement__Group_1__0__Impl rule__ShapeElement__Group_1__1 ;
    public final void rule__ShapeElement__Group_1__0() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:5349:1: ( rule__ShapeElement__Group_1__0__Impl rule__ShapeElement__Group_1__1 )
            // InternalElkGraphJson.g:5350:2: rule__ShapeElement__Group_1__0__Impl rule__ShapeElement__Group_1__1
            {
            pushFollow(FOLLOW_8);
            rule__ShapeElement__Group_1__0__Impl();

            state._fsp--;

            pushFollow(FOLLOW_2);
            rule__ShapeElement__Group_1__1();

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
    // $ANTLR end "rule__ShapeElement__Group_1__0"


    // $ANTLR start "rule__ShapeElement__Group_1__0__Impl"
    // InternalElkGraphJson.g:5357:1: rule__ShapeElement__Group_1__0__Impl : ( ruleKeyY ) ;
    public final void rule__ShapeElement__Group_1__0__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:5361:1: ( ( ruleKeyY ) )
            // InternalElkGraphJson.g:5362:1: ( ruleKeyY )
            {
            // InternalElkGraphJson.g:5362:1: ( ruleKeyY )
            // InternalElkGraphJson.g:5363:2: ruleKeyY
            {
             before(grammarAccess.getShapeElementAccess().getKeyYParserRuleCall_1_0()); 
            pushFollow(FOLLOW_2);
            ruleKeyY();

            state._fsp--;

             after(grammarAccess.getShapeElementAccess().getKeyYParserRuleCall_1_0()); 

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
    // $ANTLR end "rule__ShapeElement__Group_1__0__Impl"


    // $ANTLR start "rule__ShapeElement__Group_1__1"
    // InternalElkGraphJson.g:5372:1: rule__ShapeElement__Group_1__1 : rule__ShapeElement__Group_1__1__Impl rule__ShapeElement__Group_1__2 ;
    public final void rule__ShapeElement__Group_1__1() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:5376:1: ( rule__ShapeElement__Group_1__1__Impl rule__ShapeElement__Group_1__2 )
            // InternalElkGraphJson.g:5377:2: rule__ShapeElement__Group_1__1__Impl rule__ShapeElement__Group_1__2
            {
            pushFollow(FOLLOW_18);
            rule__ShapeElement__Group_1__1__Impl();

            state._fsp--;

            pushFollow(FOLLOW_2);
            rule__ShapeElement__Group_1__2();

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
    // $ANTLR end "rule__ShapeElement__Group_1__1"


    // $ANTLR start "rule__ShapeElement__Group_1__1__Impl"
    // InternalElkGraphJson.g:5384:1: rule__ShapeElement__Group_1__1__Impl : ( ':' ) ;
    public final void rule__ShapeElement__Group_1__1__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:5388:1: ( ( ':' ) )
            // InternalElkGraphJson.g:5389:1: ( ':' )
            {
            // InternalElkGraphJson.g:5389:1: ( ':' )
            // InternalElkGraphJson.g:5390:2: ':'
            {
             before(grammarAccess.getShapeElementAccess().getColonKeyword_1_1()); 
            match(input,61,FOLLOW_2); 
             after(grammarAccess.getShapeElementAccess().getColonKeyword_1_1()); 

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
    // $ANTLR end "rule__ShapeElement__Group_1__1__Impl"


    // $ANTLR start "rule__ShapeElement__Group_1__2"
    // InternalElkGraphJson.g:5399:1: rule__ShapeElement__Group_1__2 : rule__ShapeElement__Group_1__2__Impl ;
    public final void rule__ShapeElement__Group_1__2() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:5403:1: ( rule__ShapeElement__Group_1__2__Impl )
            // InternalElkGraphJson.g:5404:2: rule__ShapeElement__Group_1__2__Impl
            {
            pushFollow(FOLLOW_2);
            rule__ShapeElement__Group_1__2__Impl();

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
    // $ANTLR end "rule__ShapeElement__Group_1__2"


    // $ANTLR start "rule__ShapeElement__Group_1__2__Impl"
    // InternalElkGraphJson.g:5410:1: rule__ShapeElement__Group_1__2__Impl : ( ( rule__ShapeElement__YAssignment_1_2 ) ) ;
    public final void rule__ShapeElement__Group_1__2__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:5414:1: ( ( ( rule__ShapeElement__YAssignment_1_2 ) ) )
            // InternalElkGraphJson.g:5415:1: ( ( rule__ShapeElement__YAssignment_1_2 ) )
            {
            // InternalElkGraphJson.g:5415:1: ( ( rule__ShapeElement__YAssignment_1_2 ) )
            // InternalElkGraphJson.g:5416:2: ( rule__ShapeElement__YAssignment_1_2 )
            {
             before(grammarAccess.getShapeElementAccess().getYAssignment_1_2()); 
            // InternalElkGraphJson.g:5417:2: ( rule__ShapeElement__YAssignment_1_2 )
            // InternalElkGraphJson.g:5417:3: rule__ShapeElement__YAssignment_1_2
            {
            pushFollow(FOLLOW_2);
            rule__ShapeElement__YAssignment_1_2();

            state._fsp--;


            }

             after(grammarAccess.getShapeElementAccess().getYAssignment_1_2()); 

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
    // $ANTLR end "rule__ShapeElement__Group_1__2__Impl"


    // $ANTLR start "rule__ShapeElement__Group_2__0"
    // InternalElkGraphJson.g:5426:1: rule__ShapeElement__Group_2__0 : rule__ShapeElement__Group_2__0__Impl rule__ShapeElement__Group_2__1 ;
    public final void rule__ShapeElement__Group_2__0() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:5430:1: ( rule__ShapeElement__Group_2__0__Impl rule__ShapeElement__Group_2__1 )
            // InternalElkGraphJson.g:5431:2: rule__ShapeElement__Group_2__0__Impl rule__ShapeElement__Group_2__1
            {
            pushFollow(FOLLOW_8);
            rule__ShapeElement__Group_2__0__Impl();

            state._fsp--;

            pushFollow(FOLLOW_2);
            rule__ShapeElement__Group_2__1();

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
    // $ANTLR end "rule__ShapeElement__Group_2__0"


    // $ANTLR start "rule__ShapeElement__Group_2__0__Impl"
    // InternalElkGraphJson.g:5438:1: rule__ShapeElement__Group_2__0__Impl : ( ruleKeyWidth ) ;
    public final void rule__ShapeElement__Group_2__0__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:5442:1: ( ( ruleKeyWidth ) )
            // InternalElkGraphJson.g:5443:1: ( ruleKeyWidth )
            {
            // InternalElkGraphJson.g:5443:1: ( ruleKeyWidth )
            // InternalElkGraphJson.g:5444:2: ruleKeyWidth
            {
             before(grammarAccess.getShapeElementAccess().getKeyWidthParserRuleCall_2_0()); 
            pushFollow(FOLLOW_2);
            ruleKeyWidth();

            state._fsp--;

             after(grammarAccess.getShapeElementAccess().getKeyWidthParserRuleCall_2_0()); 

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
    // $ANTLR end "rule__ShapeElement__Group_2__0__Impl"


    // $ANTLR start "rule__ShapeElement__Group_2__1"
    // InternalElkGraphJson.g:5453:1: rule__ShapeElement__Group_2__1 : rule__ShapeElement__Group_2__1__Impl rule__ShapeElement__Group_2__2 ;
    public final void rule__ShapeElement__Group_2__1() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:5457:1: ( rule__ShapeElement__Group_2__1__Impl rule__ShapeElement__Group_2__2 )
            // InternalElkGraphJson.g:5458:2: rule__ShapeElement__Group_2__1__Impl rule__ShapeElement__Group_2__2
            {
            pushFollow(FOLLOW_18);
            rule__ShapeElement__Group_2__1__Impl();

            state._fsp--;

            pushFollow(FOLLOW_2);
            rule__ShapeElement__Group_2__2();

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
    // $ANTLR end "rule__ShapeElement__Group_2__1"


    // $ANTLR start "rule__ShapeElement__Group_2__1__Impl"
    // InternalElkGraphJson.g:5465:1: rule__ShapeElement__Group_2__1__Impl : ( ':' ) ;
    public final void rule__ShapeElement__Group_2__1__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:5469:1: ( ( ':' ) )
            // InternalElkGraphJson.g:5470:1: ( ':' )
            {
            // InternalElkGraphJson.g:5470:1: ( ':' )
            // InternalElkGraphJson.g:5471:2: ':'
            {
             before(grammarAccess.getShapeElementAccess().getColonKeyword_2_1()); 
            match(input,61,FOLLOW_2); 
             after(grammarAccess.getShapeElementAccess().getColonKeyword_2_1()); 

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
    // $ANTLR end "rule__ShapeElement__Group_2__1__Impl"


    // $ANTLR start "rule__ShapeElement__Group_2__2"
    // InternalElkGraphJson.g:5480:1: rule__ShapeElement__Group_2__2 : rule__ShapeElement__Group_2__2__Impl ;
    public final void rule__ShapeElement__Group_2__2() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:5484:1: ( rule__ShapeElement__Group_2__2__Impl )
            // InternalElkGraphJson.g:5485:2: rule__ShapeElement__Group_2__2__Impl
            {
            pushFollow(FOLLOW_2);
            rule__ShapeElement__Group_2__2__Impl();

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
    // $ANTLR end "rule__ShapeElement__Group_2__2"


    // $ANTLR start "rule__ShapeElement__Group_2__2__Impl"
    // InternalElkGraphJson.g:5491:1: rule__ShapeElement__Group_2__2__Impl : ( ( rule__ShapeElement__WidthAssignment_2_2 ) ) ;
    public final void rule__ShapeElement__Group_2__2__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:5495:1: ( ( ( rule__ShapeElement__WidthAssignment_2_2 ) ) )
            // InternalElkGraphJson.g:5496:1: ( ( rule__ShapeElement__WidthAssignment_2_2 ) )
            {
            // InternalElkGraphJson.g:5496:1: ( ( rule__ShapeElement__WidthAssignment_2_2 ) )
            // InternalElkGraphJson.g:5497:2: ( rule__ShapeElement__WidthAssignment_2_2 )
            {
             before(grammarAccess.getShapeElementAccess().getWidthAssignment_2_2()); 
            // InternalElkGraphJson.g:5498:2: ( rule__ShapeElement__WidthAssignment_2_2 )
            // InternalElkGraphJson.g:5498:3: rule__ShapeElement__WidthAssignment_2_2
            {
            pushFollow(FOLLOW_2);
            rule__ShapeElement__WidthAssignment_2_2();

            state._fsp--;


            }

             after(grammarAccess.getShapeElementAccess().getWidthAssignment_2_2()); 

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
    // $ANTLR end "rule__ShapeElement__Group_2__2__Impl"


    // $ANTLR start "rule__ShapeElement__Group_3__0"
    // InternalElkGraphJson.g:5507:1: rule__ShapeElement__Group_3__0 : rule__ShapeElement__Group_3__0__Impl rule__ShapeElement__Group_3__1 ;
    public final void rule__ShapeElement__Group_3__0() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:5511:1: ( rule__ShapeElement__Group_3__0__Impl rule__ShapeElement__Group_3__1 )
            // InternalElkGraphJson.g:5512:2: rule__ShapeElement__Group_3__0__Impl rule__ShapeElement__Group_3__1
            {
            pushFollow(FOLLOW_8);
            rule__ShapeElement__Group_3__0__Impl();

            state._fsp--;

            pushFollow(FOLLOW_2);
            rule__ShapeElement__Group_3__1();

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
    // $ANTLR end "rule__ShapeElement__Group_3__0"


    // $ANTLR start "rule__ShapeElement__Group_3__0__Impl"
    // InternalElkGraphJson.g:5519:1: rule__ShapeElement__Group_3__0__Impl : ( ruleKeyHeight ) ;
    public final void rule__ShapeElement__Group_3__0__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:5523:1: ( ( ruleKeyHeight ) )
            // InternalElkGraphJson.g:5524:1: ( ruleKeyHeight )
            {
            // InternalElkGraphJson.g:5524:1: ( ruleKeyHeight )
            // InternalElkGraphJson.g:5525:2: ruleKeyHeight
            {
             before(grammarAccess.getShapeElementAccess().getKeyHeightParserRuleCall_3_0()); 
            pushFollow(FOLLOW_2);
            ruleKeyHeight();

            state._fsp--;

             after(grammarAccess.getShapeElementAccess().getKeyHeightParserRuleCall_3_0()); 

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
    // $ANTLR end "rule__ShapeElement__Group_3__0__Impl"


    // $ANTLR start "rule__ShapeElement__Group_3__1"
    // InternalElkGraphJson.g:5534:1: rule__ShapeElement__Group_3__1 : rule__ShapeElement__Group_3__1__Impl rule__ShapeElement__Group_3__2 ;
    public final void rule__ShapeElement__Group_3__1() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:5538:1: ( rule__ShapeElement__Group_3__1__Impl rule__ShapeElement__Group_3__2 )
            // InternalElkGraphJson.g:5539:2: rule__ShapeElement__Group_3__1__Impl rule__ShapeElement__Group_3__2
            {
            pushFollow(FOLLOW_18);
            rule__ShapeElement__Group_3__1__Impl();

            state._fsp--;

            pushFollow(FOLLOW_2);
            rule__ShapeElement__Group_3__2();

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
    // $ANTLR end "rule__ShapeElement__Group_3__1"


    // $ANTLR start "rule__ShapeElement__Group_3__1__Impl"
    // InternalElkGraphJson.g:5546:1: rule__ShapeElement__Group_3__1__Impl : ( ':' ) ;
    public final void rule__ShapeElement__Group_3__1__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:5550:1: ( ( ':' ) )
            // InternalElkGraphJson.g:5551:1: ( ':' )
            {
            // InternalElkGraphJson.g:5551:1: ( ':' )
            // InternalElkGraphJson.g:5552:2: ':'
            {
             before(grammarAccess.getShapeElementAccess().getColonKeyword_3_1()); 
            match(input,61,FOLLOW_2); 
             after(grammarAccess.getShapeElementAccess().getColonKeyword_3_1()); 

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
    // $ANTLR end "rule__ShapeElement__Group_3__1__Impl"


    // $ANTLR start "rule__ShapeElement__Group_3__2"
    // InternalElkGraphJson.g:5561:1: rule__ShapeElement__Group_3__2 : rule__ShapeElement__Group_3__2__Impl ;
    public final void rule__ShapeElement__Group_3__2() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:5565:1: ( rule__ShapeElement__Group_3__2__Impl )
            // InternalElkGraphJson.g:5566:2: rule__ShapeElement__Group_3__2__Impl
            {
            pushFollow(FOLLOW_2);
            rule__ShapeElement__Group_3__2__Impl();

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
    // $ANTLR end "rule__ShapeElement__Group_3__2"


    // $ANTLR start "rule__ShapeElement__Group_3__2__Impl"
    // InternalElkGraphJson.g:5572:1: rule__ShapeElement__Group_3__2__Impl : ( ( rule__ShapeElement__HeightAssignment_3_2 ) ) ;
    public final void rule__ShapeElement__Group_3__2__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:5576:1: ( ( ( rule__ShapeElement__HeightAssignment_3_2 ) ) )
            // InternalElkGraphJson.g:5577:1: ( ( rule__ShapeElement__HeightAssignment_3_2 ) )
            {
            // InternalElkGraphJson.g:5577:1: ( ( rule__ShapeElement__HeightAssignment_3_2 ) )
            // InternalElkGraphJson.g:5578:2: ( rule__ShapeElement__HeightAssignment_3_2 )
            {
             before(grammarAccess.getShapeElementAccess().getHeightAssignment_3_2()); 
            // InternalElkGraphJson.g:5579:2: ( rule__ShapeElement__HeightAssignment_3_2 )
            // InternalElkGraphJson.g:5579:3: rule__ShapeElement__HeightAssignment_3_2
            {
            pushFollow(FOLLOW_2);
            rule__ShapeElement__HeightAssignment_3_2();

            state._fsp--;


            }

             after(grammarAccess.getShapeElementAccess().getHeightAssignment_3_2()); 

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
    // $ANTLR end "rule__ShapeElement__Group_3__2__Impl"


    // $ANTLR start "rule__Property__Group__0"
    // InternalElkGraphJson.g:5588:1: rule__Property__Group__0 : rule__Property__Group__0__Impl rule__Property__Group__1 ;
    public final void rule__Property__Group__0() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:5592:1: ( rule__Property__Group__0__Impl rule__Property__Group__1 )
            // InternalElkGraphJson.g:5593:2: rule__Property__Group__0__Impl rule__Property__Group__1
            {
            pushFollow(FOLLOW_8);
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
    // InternalElkGraphJson.g:5600:1: rule__Property__Group__0__Impl : ( ( rule__Property__KeyAssignment_0 ) ) ;
    public final void rule__Property__Group__0__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:5604:1: ( ( ( rule__Property__KeyAssignment_0 ) ) )
            // InternalElkGraphJson.g:5605:1: ( ( rule__Property__KeyAssignment_0 ) )
            {
            // InternalElkGraphJson.g:5605:1: ( ( rule__Property__KeyAssignment_0 ) )
            // InternalElkGraphJson.g:5606:2: ( rule__Property__KeyAssignment_0 )
            {
             before(grammarAccess.getPropertyAccess().getKeyAssignment_0()); 
            // InternalElkGraphJson.g:5607:2: ( rule__Property__KeyAssignment_0 )
            // InternalElkGraphJson.g:5607:3: rule__Property__KeyAssignment_0
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
    // InternalElkGraphJson.g:5615:1: rule__Property__Group__1 : rule__Property__Group__1__Impl rule__Property__Group__2 ;
    public final void rule__Property__Group__1() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:5619:1: ( rule__Property__Group__1__Impl rule__Property__Group__2 )
            // InternalElkGraphJson.g:5620:2: rule__Property__Group__1__Impl rule__Property__Group__2
            {
            pushFollow(FOLLOW_19);
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
    // InternalElkGraphJson.g:5627:1: rule__Property__Group__1__Impl : ( ':' ) ;
    public final void rule__Property__Group__1__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:5631:1: ( ( ':' ) )
            // InternalElkGraphJson.g:5632:1: ( ':' )
            {
            // InternalElkGraphJson.g:5632:1: ( ':' )
            // InternalElkGraphJson.g:5633:2: ':'
            {
             before(grammarAccess.getPropertyAccess().getColonKeyword_1()); 
            match(input,61,FOLLOW_2); 
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
    // InternalElkGraphJson.g:5642:1: rule__Property__Group__2 : rule__Property__Group__2__Impl ;
    public final void rule__Property__Group__2() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:5646:1: ( rule__Property__Group__2__Impl )
            // InternalElkGraphJson.g:5647:2: rule__Property__Group__2__Impl
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
    // InternalElkGraphJson.g:5653:1: rule__Property__Group__2__Impl : ( ( rule__Property__Alternatives_2 ) ) ;
    public final void rule__Property__Group__2__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:5657:1: ( ( ( rule__Property__Alternatives_2 ) ) )
            // InternalElkGraphJson.g:5658:1: ( ( rule__Property__Alternatives_2 ) )
            {
            // InternalElkGraphJson.g:5658:1: ( ( rule__Property__Alternatives_2 ) )
            // InternalElkGraphJson.g:5659:2: ( rule__Property__Alternatives_2 )
            {
             before(grammarAccess.getPropertyAccess().getAlternatives_2()); 
            // InternalElkGraphJson.g:5660:2: ( rule__Property__Alternatives_2 )
            // InternalElkGraphJson.g:5660:3: rule__Property__Alternatives_2
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


    // $ANTLR start "rule__JsonObject__Group__0"
    // InternalElkGraphJson.g:5669:1: rule__JsonObject__Group__0 : rule__JsonObject__Group__0__Impl rule__JsonObject__Group__1 ;
    public final void rule__JsonObject__Group__0() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:5673:1: ( rule__JsonObject__Group__0__Impl rule__JsonObject__Group__1 )
            // InternalElkGraphJson.g:5674:2: rule__JsonObject__Group__0__Impl rule__JsonObject__Group__1
            {
            pushFollow(FOLLOW_4);
            rule__JsonObject__Group__0__Impl();

            state._fsp--;

            pushFollow(FOLLOW_2);
            rule__JsonObject__Group__1();

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
    // $ANTLR end "rule__JsonObject__Group__0"


    // $ANTLR start "rule__JsonObject__Group__0__Impl"
    // InternalElkGraphJson.g:5681:1: rule__JsonObject__Group__0__Impl : ( '{' ) ;
    public final void rule__JsonObject__Group__0__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:5685:1: ( ( '{' ) )
            // InternalElkGraphJson.g:5686:1: ( '{' )
            {
            // InternalElkGraphJson.g:5686:1: ( '{' )
            // InternalElkGraphJson.g:5687:2: '{'
            {
             before(grammarAccess.getJsonObjectAccess().getLeftCurlyBracketKeyword_0()); 
            match(input,58,FOLLOW_2); 
             after(grammarAccess.getJsonObjectAccess().getLeftCurlyBracketKeyword_0()); 

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
    // $ANTLR end "rule__JsonObject__Group__0__Impl"


    // $ANTLR start "rule__JsonObject__Group__1"
    // InternalElkGraphJson.g:5696:1: rule__JsonObject__Group__1 : rule__JsonObject__Group__1__Impl rule__JsonObject__Group__2 ;
    public final void rule__JsonObject__Group__1() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:5700:1: ( rule__JsonObject__Group__1__Impl rule__JsonObject__Group__2 )
            // InternalElkGraphJson.g:5701:2: rule__JsonObject__Group__1__Impl rule__JsonObject__Group__2
            {
            pushFollow(FOLLOW_4);
            rule__JsonObject__Group__1__Impl();

            state._fsp--;

            pushFollow(FOLLOW_2);
            rule__JsonObject__Group__2();

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
    // $ANTLR end "rule__JsonObject__Group__1"


    // $ANTLR start "rule__JsonObject__Group__1__Impl"
    // InternalElkGraphJson.g:5708:1: rule__JsonObject__Group__1__Impl : ( ( rule__JsonObject__Group_1__0 )? ) ;
    public final void rule__JsonObject__Group__1__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:5712:1: ( ( ( rule__JsonObject__Group_1__0 )? ) )
            // InternalElkGraphJson.g:5713:1: ( ( rule__JsonObject__Group_1__0 )? )
            {
            // InternalElkGraphJson.g:5713:1: ( ( rule__JsonObject__Group_1__0 )? )
            // InternalElkGraphJson.g:5714:2: ( rule__JsonObject__Group_1__0 )?
            {
             before(grammarAccess.getJsonObjectAccess().getGroup_1()); 
            // InternalElkGraphJson.g:5715:2: ( rule__JsonObject__Group_1__0 )?
            int alt56=2;
            int LA56_0 = input.LA(1);

            if ( ((LA56_0>=RULE_STRING && LA56_0<=RULE_ID)) ) {
                alt56=1;
            }
            switch (alt56) {
                case 1 :
                    // InternalElkGraphJson.g:5715:3: rule__JsonObject__Group_1__0
                    {
                    pushFollow(FOLLOW_2);
                    rule__JsonObject__Group_1__0();

                    state._fsp--;


                    }
                    break;

            }

             after(grammarAccess.getJsonObjectAccess().getGroup_1()); 

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
    // $ANTLR end "rule__JsonObject__Group__1__Impl"


    // $ANTLR start "rule__JsonObject__Group__2"
    // InternalElkGraphJson.g:5723:1: rule__JsonObject__Group__2 : rule__JsonObject__Group__2__Impl rule__JsonObject__Group__3 ;
    public final void rule__JsonObject__Group__2() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:5727:1: ( rule__JsonObject__Group__2__Impl rule__JsonObject__Group__3 )
            // InternalElkGraphJson.g:5728:2: rule__JsonObject__Group__2__Impl rule__JsonObject__Group__3
            {
            pushFollow(FOLLOW_4);
            rule__JsonObject__Group__2__Impl();

            state._fsp--;

            pushFollow(FOLLOW_2);
            rule__JsonObject__Group__3();

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
    // $ANTLR end "rule__JsonObject__Group__2"


    // $ANTLR start "rule__JsonObject__Group__2__Impl"
    // InternalElkGraphJson.g:5735:1: rule__JsonObject__Group__2__Impl : ( ( ',' )? ) ;
    public final void rule__JsonObject__Group__2__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:5739:1: ( ( ( ',' )? ) )
            // InternalElkGraphJson.g:5740:1: ( ( ',' )? )
            {
            // InternalElkGraphJson.g:5740:1: ( ( ',' )? )
            // InternalElkGraphJson.g:5741:2: ( ',' )?
            {
             before(grammarAccess.getJsonObjectAccess().getCommaKeyword_2()); 
            // InternalElkGraphJson.g:5742:2: ( ',' )?
            int alt57=2;
            int LA57_0 = input.LA(1);

            if ( (LA57_0==59) ) {
                alt57=1;
            }
            switch (alt57) {
                case 1 :
                    // InternalElkGraphJson.g:5742:3: ','
                    {
                    match(input,59,FOLLOW_2); 

                    }
                    break;

            }

             after(grammarAccess.getJsonObjectAccess().getCommaKeyword_2()); 

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
    // $ANTLR end "rule__JsonObject__Group__2__Impl"


    // $ANTLR start "rule__JsonObject__Group__3"
    // InternalElkGraphJson.g:5750:1: rule__JsonObject__Group__3 : rule__JsonObject__Group__3__Impl ;
    public final void rule__JsonObject__Group__3() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:5754:1: ( rule__JsonObject__Group__3__Impl )
            // InternalElkGraphJson.g:5755:2: rule__JsonObject__Group__3__Impl
            {
            pushFollow(FOLLOW_2);
            rule__JsonObject__Group__3__Impl();

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
    // $ANTLR end "rule__JsonObject__Group__3"


    // $ANTLR start "rule__JsonObject__Group__3__Impl"
    // InternalElkGraphJson.g:5761:1: rule__JsonObject__Group__3__Impl : ( '}' ) ;
    public final void rule__JsonObject__Group__3__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:5765:1: ( ( '}' ) )
            // InternalElkGraphJson.g:5766:1: ( '}' )
            {
            // InternalElkGraphJson.g:5766:1: ( '}' )
            // InternalElkGraphJson.g:5767:2: '}'
            {
             before(grammarAccess.getJsonObjectAccess().getRightCurlyBracketKeyword_3()); 
            match(input,60,FOLLOW_2); 
             after(grammarAccess.getJsonObjectAccess().getRightCurlyBracketKeyword_3()); 

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
    // $ANTLR end "rule__JsonObject__Group__3__Impl"


    // $ANTLR start "rule__JsonObject__Group_1__0"
    // InternalElkGraphJson.g:5777:1: rule__JsonObject__Group_1__0 : rule__JsonObject__Group_1__0__Impl rule__JsonObject__Group_1__1 ;
    public final void rule__JsonObject__Group_1__0() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:5781:1: ( rule__JsonObject__Group_1__0__Impl rule__JsonObject__Group_1__1 )
            // InternalElkGraphJson.g:5782:2: rule__JsonObject__Group_1__0__Impl rule__JsonObject__Group_1__1
            {
            pushFollow(FOLLOW_5);
            rule__JsonObject__Group_1__0__Impl();

            state._fsp--;

            pushFollow(FOLLOW_2);
            rule__JsonObject__Group_1__1();

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
    // $ANTLR end "rule__JsonObject__Group_1__0"


    // $ANTLR start "rule__JsonObject__Group_1__0__Impl"
    // InternalElkGraphJson.g:5789:1: rule__JsonObject__Group_1__0__Impl : ( ruleJsonMember ) ;
    public final void rule__JsonObject__Group_1__0__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:5793:1: ( ( ruleJsonMember ) )
            // InternalElkGraphJson.g:5794:1: ( ruleJsonMember )
            {
            // InternalElkGraphJson.g:5794:1: ( ruleJsonMember )
            // InternalElkGraphJson.g:5795:2: ruleJsonMember
            {
             before(grammarAccess.getJsonObjectAccess().getJsonMemberParserRuleCall_1_0()); 
            pushFollow(FOLLOW_2);
            ruleJsonMember();

            state._fsp--;

             after(grammarAccess.getJsonObjectAccess().getJsonMemberParserRuleCall_1_0()); 

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
    // $ANTLR end "rule__JsonObject__Group_1__0__Impl"


    // $ANTLR start "rule__JsonObject__Group_1__1"
    // InternalElkGraphJson.g:5804:1: rule__JsonObject__Group_1__1 : rule__JsonObject__Group_1__1__Impl ;
    public final void rule__JsonObject__Group_1__1() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:5808:1: ( rule__JsonObject__Group_1__1__Impl )
            // InternalElkGraphJson.g:5809:2: rule__JsonObject__Group_1__1__Impl
            {
            pushFollow(FOLLOW_2);
            rule__JsonObject__Group_1__1__Impl();

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
    // $ANTLR end "rule__JsonObject__Group_1__1"


    // $ANTLR start "rule__JsonObject__Group_1__1__Impl"
    // InternalElkGraphJson.g:5815:1: rule__JsonObject__Group_1__1__Impl : ( ( rule__JsonObject__Group_1_1__0 )* ) ;
    public final void rule__JsonObject__Group_1__1__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:5819:1: ( ( ( rule__JsonObject__Group_1_1__0 )* ) )
            // InternalElkGraphJson.g:5820:1: ( ( rule__JsonObject__Group_1_1__0 )* )
            {
            // InternalElkGraphJson.g:5820:1: ( ( rule__JsonObject__Group_1_1__0 )* )
            // InternalElkGraphJson.g:5821:2: ( rule__JsonObject__Group_1_1__0 )*
            {
             before(grammarAccess.getJsonObjectAccess().getGroup_1_1()); 
            // InternalElkGraphJson.g:5822:2: ( rule__JsonObject__Group_1_1__0 )*
            loop58:
            do {
                int alt58=2;
                int LA58_0 = input.LA(1);

                if ( (LA58_0==59) ) {
                    int LA58_1 = input.LA(2);

                    if ( ((LA58_1>=RULE_STRING && LA58_1<=RULE_ID)) ) {
                        alt58=1;
                    }


                }


                switch (alt58) {
            	case 1 :
            	    // InternalElkGraphJson.g:5822:3: rule__JsonObject__Group_1_1__0
            	    {
            	    pushFollow(FOLLOW_6);
            	    rule__JsonObject__Group_1_1__0();

            	    state._fsp--;


            	    }
            	    break;

            	default :
            	    break loop58;
                }
            } while (true);

             after(grammarAccess.getJsonObjectAccess().getGroup_1_1()); 

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
    // $ANTLR end "rule__JsonObject__Group_1__1__Impl"


    // $ANTLR start "rule__JsonObject__Group_1_1__0"
    // InternalElkGraphJson.g:5831:1: rule__JsonObject__Group_1_1__0 : rule__JsonObject__Group_1_1__0__Impl rule__JsonObject__Group_1_1__1 ;
    public final void rule__JsonObject__Group_1_1__0() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:5835:1: ( rule__JsonObject__Group_1_1__0__Impl rule__JsonObject__Group_1_1__1 )
            // InternalElkGraphJson.g:5836:2: rule__JsonObject__Group_1_1__0__Impl rule__JsonObject__Group_1_1__1
            {
            pushFollow(FOLLOW_7);
            rule__JsonObject__Group_1_1__0__Impl();

            state._fsp--;

            pushFollow(FOLLOW_2);
            rule__JsonObject__Group_1_1__1();

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
    // $ANTLR end "rule__JsonObject__Group_1_1__0"


    // $ANTLR start "rule__JsonObject__Group_1_1__0__Impl"
    // InternalElkGraphJson.g:5843:1: rule__JsonObject__Group_1_1__0__Impl : ( ',' ) ;
    public final void rule__JsonObject__Group_1_1__0__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:5847:1: ( ( ',' ) )
            // InternalElkGraphJson.g:5848:1: ( ',' )
            {
            // InternalElkGraphJson.g:5848:1: ( ',' )
            // InternalElkGraphJson.g:5849:2: ','
            {
             before(grammarAccess.getJsonObjectAccess().getCommaKeyword_1_1_0()); 
            match(input,59,FOLLOW_2); 
             after(grammarAccess.getJsonObjectAccess().getCommaKeyword_1_1_0()); 

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
    // $ANTLR end "rule__JsonObject__Group_1_1__0__Impl"


    // $ANTLR start "rule__JsonObject__Group_1_1__1"
    // InternalElkGraphJson.g:5858:1: rule__JsonObject__Group_1_1__1 : rule__JsonObject__Group_1_1__1__Impl ;
    public final void rule__JsonObject__Group_1_1__1() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:5862:1: ( rule__JsonObject__Group_1_1__1__Impl )
            // InternalElkGraphJson.g:5863:2: rule__JsonObject__Group_1_1__1__Impl
            {
            pushFollow(FOLLOW_2);
            rule__JsonObject__Group_1_1__1__Impl();

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
    // $ANTLR end "rule__JsonObject__Group_1_1__1"


    // $ANTLR start "rule__JsonObject__Group_1_1__1__Impl"
    // InternalElkGraphJson.g:5869:1: rule__JsonObject__Group_1_1__1__Impl : ( ruleJsonMember ) ;
    public final void rule__JsonObject__Group_1_1__1__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:5873:1: ( ( ruleJsonMember ) )
            // InternalElkGraphJson.g:5874:1: ( ruleJsonMember )
            {
            // InternalElkGraphJson.g:5874:1: ( ruleJsonMember )
            // InternalElkGraphJson.g:5875:2: ruleJsonMember
            {
             before(grammarAccess.getJsonObjectAccess().getJsonMemberParserRuleCall_1_1_1()); 
            pushFollow(FOLLOW_2);
            ruleJsonMember();

            state._fsp--;

             after(grammarAccess.getJsonObjectAccess().getJsonMemberParserRuleCall_1_1_1()); 

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
    // $ANTLR end "rule__JsonObject__Group_1_1__1__Impl"


    // $ANTLR start "rule__JsonArray__Group__0"
    // InternalElkGraphJson.g:5885:1: rule__JsonArray__Group__0 : rule__JsonArray__Group__0__Impl rule__JsonArray__Group__1 ;
    public final void rule__JsonArray__Group__0() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:5889:1: ( rule__JsonArray__Group__0__Impl rule__JsonArray__Group__1 )
            // InternalElkGraphJson.g:5890:2: rule__JsonArray__Group__0__Impl rule__JsonArray__Group__1
            {
            pushFollow(FOLLOW_20);
            rule__JsonArray__Group__0__Impl();

            state._fsp--;

            pushFollow(FOLLOW_2);
            rule__JsonArray__Group__1();

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
    // $ANTLR end "rule__JsonArray__Group__0"


    // $ANTLR start "rule__JsonArray__Group__0__Impl"
    // InternalElkGraphJson.g:5897:1: rule__JsonArray__Group__0__Impl : ( '[' ) ;
    public final void rule__JsonArray__Group__0__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:5901:1: ( ( '[' ) )
            // InternalElkGraphJson.g:5902:1: ( '[' )
            {
            // InternalElkGraphJson.g:5902:1: ( '[' )
            // InternalElkGraphJson.g:5903:2: '['
            {
             before(grammarAccess.getJsonArrayAccess().getLeftSquareBracketKeyword_0()); 
            match(input,62,FOLLOW_2); 
             after(grammarAccess.getJsonArrayAccess().getLeftSquareBracketKeyword_0()); 

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
    // $ANTLR end "rule__JsonArray__Group__0__Impl"


    // $ANTLR start "rule__JsonArray__Group__1"
    // InternalElkGraphJson.g:5912:1: rule__JsonArray__Group__1 : rule__JsonArray__Group__1__Impl rule__JsonArray__Group__2 ;
    public final void rule__JsonArray__Group__1() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:5916:1: ( rule__JsonArray__Group__1__Impl rule__JsonArray__Group__2 )
            // InternalElkGraphJson.g:5917:2: rule__JsonArray__Group__1__Impl rule__JsonArray__Group__2
            {
            pushFollow(FOLLOW_20);
            rule__JsonArray__Group__1__Impl();

            state._fsp--;

            pushFollow(FOLLOW_2);
            rule__JsonArray__Group__2();

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
    // $ANTLR end "rule__JsonArray__Group__1"


    // $ANTLR start "rule__JsonArray__Group__1__Impl"
    // InternalElkGraphJson.g:5924:1: rule__JsonArray__Group__1__Impl : ( ( rule__JsonArray__Group_1__0 )? ) ;
    public final void rule__JsonArray__Group__1__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:5928:1: ( ( ( rule__JsonArray__Group_1__0 )? ) )
            // InternalElkGraphJson.g:5929:1: ( ( rule__JsonArray__Group_1__0 )? )
            {
            // InternalElkGraphJson.g:5929:1: ( ( rule__JsonArray__Group_1__0 )? )
            // InternalElkGraphJson.g:5930:2: ( rule__JsonArray__Group_1__0 )?
            {
             before(grammarAccess.getJsonArrayAccess().getGroup_1()); 
            // InternalElkGraphJson.g:5931:2: ( rule__JsonArray__Group_1__0 )?
            int alt59=2;
            int LA59_0 = input.LA(1);

            if ( (LA59_0==RULE_STRING||(LA59_0>=RULE_SIGNED_INT && LA59_0<=RULE_FLOAT)||(LA59_0>=13 && LA59_0<=15)||LA59_0==58||LA59_0==62) ) {
                alt59=1;
            }
            switch (alt59) {
                case 1 :
                    // InternalElkGraphJson.g:5931:3: rule__JsonArray__Group_1__0
                    {
                    pushFollow(FOLLOW_2);
                    rule__JsonArray__Group_1__0();

                    state._fsp--;


                    }
                    break;

            }

             after(grammarAccess.getJsonArrayAccess().getGroup_1()); 

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
    // $ANTLR end "rule__JsonArray__Group__1__Impl"


    // $ANTLR start "rule__JsonArray__Group__2"
    // InternalElkGraphJson.g:5939:1: rule__JsonArray__Group__2 : rule__JsonArray__Group__2__Impl rule__JsonArray__Group__3 ;
    public final void rule__JsonArray__Group__2() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:5943:1: ( rule__JsonArray__Group__2__Impl rule__JsonArray__Group__3 )
            // InternalElkGraphJson.g:5944:2: rule__JsonArray__Group__2__Impl rule__JsonArray__Group__3
            {
            pushFollow(FOLLOW_20);
            rule__JsonArray__Group__2__Impl();

            state._fsp--;

            pushFollow(FOLLOW_2);
            rule__JsonArray__Group__3();

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
    // $ANTLR end "rule__JsonArray__Group__2"


    // $ANTLR start "rule__JsonArray__Group__2__Impl"
    // InternalElkGraphJson.g:5951:1: rule__JsonArray__Group__2__Impl : ( ( ',' )? ) ;
    public final void rule__JsonArray__Group__2__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:5955:1: ( ( ( ',' )? ) )
            // InternalElkGraphJson.g:5956:1: ( ( ',' )? )
            {
            // InternalElkGraphJson.g:5956:1: ( ( ',' )? )
            // InternalElkGraphJson.g:5957:2: ( ',' )?
            {
             before(grammarAccess.getJsonArrayAccess().getCommaKeyword_2()); 
            // InternalElkGraphJson.g:5958:2: ( ',' )?
            int alt60=2;
            int LA60_0 = input.LA(1);

            if ( (LA60_0==59) ) {
                alt60=1;
            }
            switch (alt60) {
                case 1 :
                    // InternalElkGraphJson.g:5958:3: ','
                    {
                    match(input,59,FOLLOW_2); 

                    }
                    break;

            }

             after(grammarAccess.getJsonArrayAccess().getCommaKeyword_2()); 

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
    // $ANTLR end "rule__JsonArray__Group__2__Impl"


    // $ANTLR start "rule__JsonArray__Group__3"
    // InternalElkGraphJson.g:5966:1: rule__JsonArray__Group__3 : rule__JsonArray__Group__3__Impl ;
    public final void rule__JsonArray__Group__3() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:5970:1: ( rule__JsonArray__Group__3__Impl )
            // InternalElkGraphJson.g:5971:2: rule__JsonArray__Group__3__Impl
            {
            pushFollow(FOLLOW_2);
            rule__JsonArray__Group__3__Impl();

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
    // $ANTLR end "rule__JsonArray__Group__3"


    // $ANTLR start "rule__JsonArray__Group__3__Impl"
    // InternalElkGraphJson.g:5977:1: rule__JsonArray__Group__3__Impl : ( ']' ) ;
    public final void rule__JsonArray__Group__3__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:5981:1: ( ( ']' ) )
            // InternalElkGraphJson.g:5982:1: ( ']' )
            {
            // InternalElkGraphJson.g:5982:1: ( ']' )
            // InternalElkGraphJson.g:5983:2: ']'
            {
             before(grammarAccess.getJsonArrayAccess().getRightSquareBracketKeyword_3()); 
            match(input,63,FOLLOW_2); 
             after(grammarAccess.getJsonArrayAccess().getRightSquareBracketKeyword_3()); 

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
    // $ANTLR end "rule__JsonArray__Group__3__Impl"


    // $ANTLR start "rule__JsonArray__Group_1__0"
    // InternalElkGraphJson.g:5993:1: rule__JsonArray__Group_1__0 : rule__JsonArray__Group_1__0__Impl rule__JsonArray__Group_1__1 ;
    public final void rule__JsonArray__Group_1__0() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:5997:1: ( rule__JsonArray__Group_1__0__Impl rule__JsonArray__Group_1__1 )
            // InternalElkGraphJson.g:5998:2: rule__JsonArray__Group_1__0__Impl rule__JsonArray__Group_1__1
            {
            pushFollow(FOLLOW_5);
            rule__JsonArray__Group_1__0__Impl();

            state._fsp--;

            pushFollow(FOLLOW_2);
            rule__JsonArray__Group_1__1();

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
    // $ANTLR end "rule__JsonArray__Group_1__0"


    // $ANTLR start "rule__JsonArray__Group_1__0__Impl"
    // InternalElkGraphJson.g:6005:1: rule__JsonArray__Group_1__0__Impl : ( ruleJsonValue ) ;
    public final void rule__JsonArray__Group_1__0__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:6009:1: ( ( ruleJsonValue ) )
            // InternalElkGraphJson.g:6010:1: ( ruleJsonValue )
            {
            // InternalElkGraphJson.g:6010:1: ( ruleJsonValue )
            // InternalElkGraphJson.g:6011:2: ruleJsonValue
            {
             before(grammarAccess.getJsonArrayAccess().getJsonValueParserRuleCall_1_0()); 
            pushFollow(FOLLOW_2);
            ruleJsonValue();

            state._fsp--;

             after(grammarAccess.getJsonArrayAccess().getJsonValueParserRuleCall_1_0()); 

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
    // $ANTLR end "rule__JsonArray__Group_1__0__Impl"


    // $ANTLR start "rule__JsonArray__Group_1__1"
    // InternalElkGraphJson.g:6020:1: rule__JsonArray__Group_1__1 : rule__JsonArray__Group_1__1__Impl ;
    public final void rule__JsonArray__Group_1__1() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:6024:1: ( rule__JsonArray__Group_1__1__Impl )
            // InternalElkGraphJson.g:6025:2: rule__JsonArray__Group_1__1__Impl
            {
            pushFollow(FOLLOW_2);
            rule__JsonArray__Group_1__1__Impl();

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
    // $ANTLR end "rule__JsonArray__Group_1__1"


    // $ANTLR start "rule__JsonArray__Group_1__1__Impl"
    // InternalElkGraphJson.g:6031:1: rule__JsonArray__Group_1__1__Impl : ( ( rule__JsonArray__Group_1_1__0 )* ) ;
    public final void rule__JsonArray__Group_1__1__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:6035:1: ( ( ( rule__JsonArray__Group_1_1__0 )* ) )
            // InternalElkGraphJson.g:6036:1: ( ( rule__JsonArray__Group_1_1__0 )* )
            {
            // InternalElkGraphJson.g:6036:1: ( ( rule__JsonArray__Group_1_1__0 )* )
            // InternalElkGraphJson.g:6037:2: ( rule__JsonArray__Group_1_1__0 )*
            {
             before(grammarAccess.getJsonArrayAccess().getGroup_1_1()); 
            // InternalElkGraphJson.g:6038:2: ( rule__JsonArray__Group_1_1__0 )*
            loop61:
            do {
                int alt61=2;
                int LA61_0 = input.LA(1);

                if ( (LA61_0==59) ) {
                    int LA61_1 = input.LA(2);

                    if ( (LA61_1==RULE_STRING||(LA61_1>=RULE_SIGNED_INT && LA61_1<=RULE_FLOAT)||(LA61_1>=13 && LA61_1<=15)||LA61_1==58||LA61_1==62) ) {
                        alt61=1;
                    }


                }


                switch (alt61) {
            	case 1 :
            	    // InternalElkGraphJson.g:6038:3: rule__JsonArray__Group_1_1__0
            	    {
            	    pushFollow(FOLLOW_6);
            	    rule__JsonArray__Group_1_1__0();

            	    state._fsp--;


            	    }
            	    break;

            	default :
            	    break loop61;
                }
            } while (true);

             after(grammarAccess.getJsonArrayAccess().getGroup_1_1()); 

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
    // $ANTLR end "rule__JsonArray__Group_1__1__Impl"


    // $ANTLR start "rule__JsonArray__Group_1_1__0"
    // InternalElkGraphJson.g:6047:1: rule__JsonArray__Group_1_1__0 : rule__JsonArray__Group_1_1__0__Impl rule__JsonArray__Group_1_1__1 ;
    public final void rule__JsonArray__Group_1_1__0() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:6051:1: ( rule__JsonArray__Group_1_1__0__Impl rule__JsonArray__Group_1_1__1 )
            // InternalElkGraphJson.g:6052:2: rule__JsonArray__Group_1_1__0__Impl rule__JsonArray__Group_1_1__1
            {
            pushFollow(FOLLOW_21);
            rule__JsonArray__Group_1_1__0__Impl();

            state._fsp--;

            pushFollow(FOLLOW_2);
            rule__JsonArray__Group_1_1__1();

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
    // $ANTLR end "rule__JsonArray__Group_1_1__0"


    // $ANTLR start "rule__JsonArray__Group_1_1__0__Impl"
    // InternalElkGraphJson.g:6059:1: rule__JsonArray__Group_1_1__0__Impl : ( ',' ) ;
    public final void rule__JsonArray__Group_1_1__0__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:6063:1: ( ( ',' ) )
            // InternalElkGraphJson.g:6064:1: ( ',' )
            {
            // InternalElkGraphJson.g:6064:1: ( ',' )
            // InternalElkGraphJson.g:6065:2: ','
            {
             before(grammarAccess.getJsonArrayAccess().getCommaKeyword_1_1_0()); 
            match(input,59,FOLLOW_2); 
             after(grammarAccess.getJsonArrayAccess().getCommaKeyword_1_1_0()); 

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
    // $ANTLR end "rule__JsonArray__Group_1_1__0__Impl"


    // $ANTLR start "rule__JsonArray__Group_1_1__1"
    // InternalElkGraphJson.g:6074:1: rule__JsonArray__Group_1_1__1 : rule__JsonArray__Group_1_1__1__Impl ;
    public final void rule__JsonArray__Group_1_1__1() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:6078:1: ( rule__JsonArray__Group_1_1__1__Impl )
            // InternalElkGraphJson.g:6079:2: rule__JsonArray__Group_1_1__1__Impl
            {
            pushFollow(FOLLOW_2);
            rule__JsonArray__Group_1_1__1__Impl();

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
    // $ANTLR end "rule__JsonArray__Group_1_1__1"


    // $ANTLR start "rule__JsonArray__Group_1_1__1__Impl"
    // InternalElkGraphJson.g:6085:1: rule__JsonArray__Group_1_1__1__Impl : ( ruleJsonValue ) ;
    public final void rule__JsonArray__Group_1_1__1__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:6089:1: ( ( ruleJsonValue ) )
            // InternalElkGraphJson.g:6090:1: ( ruleJsonValue )
            {
            // InternalElkGraphJson.g:6090:1: ( ruleJsonValue )
            // InternalElkGraphJson.g:6091:2: ruleJsonValue
            {
             before(grammarAccess.getJsonArrayAccess().getJsonValueParserRuleCall_1_1_1()); 
            pushFollow(FOLLOW_2);
            ruleJsonValue();

            state._fsp--;

             after(grammarAccess.getJsonArrayAccess().getJsonValueParserRuleCall_1_1_1()); 

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
    // $ANTLR end "rule__JsonArray__Group_1_1__1__Impl"


    // $ANTLR start "rule__JsonMember__Group__0"
    // InternalElkGraphJson.g:6101:1: rule__JsonMember__Group__0 : rule__JsonMember__Group__0__Impl rule__JsonMember__Group__1 ;
    public final void rule__JsonMember__Group__0() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:6105:1: ( rule__JsonMember__Group__0__Impl rule__JsonMember__Group__1 )
            // InternalElkGraphJson.g:6106:2: rule__JsonMember__Group__0__Impl rule__JsonMember__Group__1
            {
            pushFollow(FOLLOW_8);
            rule__JsonMember__Group__0__Impl();

            state._fsp--;

            pushFollow(FOLLOW_2);
            rule__JsonMember__Group__1();

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
    // $ANTLR end "rule__JsonMember__Group__0"


    // $ANTLR start "rule__JsonMember__Group__0__Impl"
    // InternalElkGraphJson.g:6113:1: rule__JsonMember__Group__0__Impl : ( ( rule__JsonMember__Alternatives_0 ) ) ;
    public final void rule__JsonMember__Group__0__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:6117:1: ( ( ( rule__JsonMember__Alternatives_0 ) ) )
            // InternalElkGraphJson.g:6118:1: ( ( rule__JsonMember__Alternatives_0 ) )
            {
            // InternalElkGraphJson.g:6118:1: ( ( rule__JsonMember__Alternatives_0 ) )
            // InternalElkGraphJson.g:6119:2: ( rule__JsonMember__Alternatives_0 )
            {
             before(grammarAccess.getJsonMemberAccess().getAlternatives_0()); 
            // InternalElkGraphJson.g:6120:2: ( rule__JsonMember__Alternatives_0 )
            // InternalElkGraphJson.g:6120:3: rule__JsonMember__Alternatives_0
            {
            pushFollow(FOLLOW_2);
            rule__JsonMember__Alternatives_0();

            state._fsp--;


            }

             after(grammarAccess.getJsonMemberAccess().getAlternatives_0()); 

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
    // $ANTLR end "rule__JsonMember__Group__0__Impl"


    // $ANTLR start "rule__JsonMember__Group__1"
    // InternalElkGraphJson.g:6128:1: rule__JsonMember__Group__1 : rule__JsonMember__Group__1__Impl rule__JsonMember__Group__2 ;
    public final void rule__JsonMember__Group__1() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:6132:1: ( rule__JsonMember__Group__1__Impl rule__JsonMember__Group__2 )
            // InternalElkGraphJson.g:6133:2: rule__JsonMember__Group__1__Impl rule__JsonMember__Group__2
            {
            pushFollow(FOLLOW_21);
            rule__JsonMember__Group__1__Impl();

            state._fsp--;

            pushFollow(FOLLOW_2);
            rule__JsonMember__Group__2();

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
    // $ANTLR end "rule__JsonMember__Group__1"


    // $ANTLR start "rule__JsonMember__Group__1__Impl"
    // InternalElkGraphJson.g:6140:1: rule__JsonMember__Group__1__Impl : ( ':' ) ;
    public final void rule__JsonMember__Group__1__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:6144:1: ( ( ':' ) )
            // InternalElkGraphJson.g:6145:1: ( ':' )
            {
            // InternalElkGraphJson.g:6145:1: ( ':' )
            // InternalElkGraphJson.g:6146:2: ':'
            {
             before(grammarAccess.getJsonMemberAccess().getColonKeyword_1()); 
            match(input,61,FOLLOW_2); 
             after(grammarAccess.getJsonMemberAccess().getColonKeyword_1()); 

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
    // $ANTLR end "rule__JsonMember__Group__1__Impl"


    // $ANTLR start "rule__JsonMember__Group__2"
    // InternalElkGraphJson.g:6155:1: rule__JsonMember__Group__2 : rule__JsonMember__Group__2__Impl ;
    public final void rule__JsonMember__Group__2() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:6159:1: ( rule__JsonMember__Group__2__Impl )
            // InternalElkGraphJson.g:6160:2: rule__JsonMember__Group__2__Impl
            {
            pushFollow(FOLLOW_2);
            rule__JsonMember__Group__2__Impl();

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
    // $ANTLR end "rule__JsonMember__Group__2"


    // $ANTLR start "rule__JsonMember__Group__2__Impl"
    // InternalElkGraphJson.g:6166:1: rule__JsonMember__Group__2__Impl : ( ruleJsonValue ) ;
    public final void rule__JsonMember__Group__2__Impl() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:6170:1: ( ( ruleJsonValue ) )
            // InternalElkGraphJson.g:6171:1: ( ruleJsonValue )
            {
            // InternalElkGraphJson.g:6171:1: ( ruleJsonValue )
            // InternalElkGraphJson.g:6172:2: ruleJsonValue
            {
             before(grammarAccess.getJsonMemberAccess().getJsonValueParserRuleCall_2()); 
            pushFollow(FOLLOW_2);
            ruleJsonValue();

            state._fsp--;

             after(grammarAccess.getJsonMemberAccess().getJsonValueParserRuleCall_2()); 

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
    // $ANTLR end "rule__JsonMember__Group__2__Impl"


    // $ANTLR start "rule__LabelElement__TextAssignment_2_2"
    // InternalElkGraphJson.g:6182:1: rule__LabelElement__TextAssignment_2_2 : ( RULE_STRING ) ;
    public final void rule__LabelElement__TextAssignment_2_2() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:6186:1: ( ( RULE_STRING ) )
            // InternalElkGraphJson.g:6187:2: ( RULE_STRING )
            {
            // InternalElkGraphJson.g:6187:2: ( RULE_STRING )
            // InternalElkGraphJson.g:6188:3: RULE_STRING
            {
             before(grammarAccess.getLabelElementAccess().getTextSTRINGTerminalRuleCall_2_2_0()); 
            match(input,RULE_STRING,FOLLOW_2); 
             after(grammarAccess.getLabelElementAccess().getTextSTRINGTerminalRuleCall_2_2_0()); 

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
    // $ANTLR end "rule__LabelElement__TextAssignment_2_2"


    // $ANTLR start "rule__ElkEdgeSources__SourcesAssignment_1_0"
    // InternalElkGraphJson.g:6197:1: rule__ElkEdgeSources__SourcesAssignment_1_0 : ( ( RULE_STRING ) ) ;
    public final void rule__ElkEdgeSources__SourcesAssignment_1_0() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:6201:1: ( ( ( RULE_STRING ) ) )
            // InternalElkGraphJson.g:6202:2: ( ( RULE_STRING ) )
            {
            // InternalElkGraphJson.g:6202:2: ( ( RULE_STRING ) )
            // InternalElkGraphJson.g:6203:3: ( RULE_STRING )
            {
             before(grammarAccess.getElkEdgeSourcesAccess().getSourcesElkConnectableShapeCrossReference_1_0_0()); 
            // InternalElkGraphJson.g:6204:3: ( RULE_STRING )
            // InternalElkGraphJson.g:6205:4: RULE_STRING
            {
             before(grammarAccess.getElkEdgeSourcesAccess().getSourcesElkConnectableShapeSTRINGTerminalRuleCall_1_0_0_1()); 
            match(input,RULE_STRING,FOLLOW_2); 
             after(grammarAccess.getElkEdgeSourcesAccess().getSourcesElkConnectableShapeSTRINGTerminalRuleCall_1_0_0_1()); 

            }

             after(grammarAccess.getElkEdgeSourcesAccess().getSourcesElkConnectableShapeCrossReference_1_0_0()); 

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
    // $ANTLR end "rule__ElkEdgeSources__SourcesAssignment_1_0"


    // $ANTLR start "rule__ElkEdgeSources__SourcesAssignment_1_1_1"
    // InternalElkGraphJson.g:6216:1: rule__ElkEdgeSources__SourcesAssignment_1_1_1 : ( ( RULE_STRING ) ) ;
    public final void rule__ElkEdgeSources__SourcesAssignment_1_1_1() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:6220:1: ( ( ( RULE_STRING ) ) )
            // InternalElkGraphJson.g:6221:2: ( ( RULE_STRING ) )
            {
            // InternalElkGraphJson.g:6221:2: ( ( RULE_STRING ) )
            // InternalElkGraphJson.g:6222:3: ( RULE_STRING )
            {
             before(grammarAccess.getElkEdgeSourcesAccess().getSourcesElkConnectableShapeCrossReference_1_1_1_0()); 
            // InternalElkGraphJson.g:6223:3: ( RULE_STRING )
            // InternalElkGraphJson.g:6224:4: RULE_STRING
            {
             before(grammarAccess.getElkEdgeSourcesAccess().getSourcesElkConnectableShapeSTRINGTerminalRuleCall_1_1_1_0_1()); 
            match(input,RULE_STRING,FOLLOW_2); 
             after(grammarAccess.getElkEdgeSourcesAccess().getSourcesElkConnectableShapeSTRINGTerminalRuleCall_1_1_1_0_1()); 

            }

             after(grammarAccess.getElkEdgeSourcesAccess().getSourcesElkConnectableShapeCrossReference_1_1_1_0()); 

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
    // $ANTLR end "rule__ElkEdgeSources__SourcesAssignment_1_1_1"


    // $ANTLR start "rule__ElkEdgeTargets__TargetsAssignment_1_0"
    // InternalElkGraphJson.g:6235:1: rule__ElkEdgeTargets__TargetsAssignment_1_0 : ( ( RULE_STRING ) ) ;
    public final void rule__ElkEdgeTargets__TargetsAssignment_1_0() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:6239:1: ( ( ( RULE_STRING ) ) )
            // InternalElkGraphJson.g:6240:2: ( ( RULE_STRING ) )
            {
            // InternalElkGraphJson.g:6240:2: ( ( RULE_STRING ) )
            // InternalElkGraphJson.g:6241:3: ( RULE_STRING )
            {
             before(grammarAccess.getElkEdgeTargetsAccess().getTargetsElkConnectableShapeCrossReference_1_0_0()); 
            // InternalElkGraphJson.g:6242:3: ( RULE_STRING )
            // InternalElkGraphJson.g:6243:4: RULE_STRING
            {
             before(grammarAccess.getElkEdgeTargetsAccess().getTargetsElkConnectableShapeSTRINGTerminalRuleCall_1_0_0_1()); 
            match(input,RULE_STRING,FOLLOW_2); 
             after(grammarAccess.getElkEdgeTargetsAccess().getTargetsElkConnectableShapeSTRINGTerminalRuleCall_1_0_0_1()); 

            }

             after(grammarAccess.getElkEdgeTargetsAccess().getTargetsElkConnectableShapeCrossReference_1_0_0()); 

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
    // $ANTLR end "rule__ElkEdgeTargets__TargetsAssignment_1_0"


    // $ANTLR start "rule__ElkEdgeTargets__TargetsAssignment_1_1_1"
    // InternalElkGraphJson.g:6254:1: rule__ElkEdgeTargets__TargetsAssignment_1_1_1 : ( ( RULE_STRING ) ) ;
    public final void rule__ElkEdgeTargets__TargetsAssignment_1_1_1() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:6258:1: ( ( ( RULE_STRING ) ) )
            // InternalElkGraphJson.g:6259:2: ( ( RULE_STRING ) )
            {
            // InternalElkGraphJson.g:6259:2: ( ( RULE_STRING ) )
            // InternalElkGraphJson.g:6260:3: ( RULE_STRING )
            {
             before(grammarAccess.getElkEdgeTargetsAccess().getTargetsElkConnectableShapeCrossReference_1_1_1_0()); 
            // InternalElkGraphJson.g:6261:3: ( RULE_STRING )
            // InternalElkGraphJson.g:6262:4: RULE_STRING
            {
             before(grammarAccess.getElkEdgeTargetsAccess().getTargetsElkConnectableShapeSTRINGTerminalRuleCall_1_1_1_0_1()); 
            match(input,RULE_STRING,FOLLOW_2); 
             after(grammarAccess.getElkEdgeTargetsAccess().getTargetsElkConnectableShapeSTRINGTerminalRuleCall_1_1_1_0_1()); 

            }

             after(grammarAccess.getElkEdgeTargetsAccess().getTargetsElkConnectableShapeCrossReference_1_1_1_0()); 

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
    // $ANTLR end "rule__ElkEdgeTargets__TargetsAssignment_1_1_1"


    // $ANTLR start "rule__ElkId__IdentifierAssignment_2"
    // InternalElkGraphJson.g:6273:1: rule__ElkId__IdentifierAssignment_2 : ( RULE_STRING ) ;
    public final void rule__ElkId__IdentifierAssignment_2() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:6277:1: ( ( RULE_STRING ) )
            // InternalElkGraphJson.g:6278:2: ( RULE_STRING )
            {
            // InternalElkGraphJson.g:6278:2: ( RULE_STRING )
            // InternalElkGraphJson.g:6279:3: RULE_STRING
            {
             before(grammarAccess.getElkIdAccess().getIdentifierSTRINGTerminalRuleCall_2_0()); 
            match(input,RULE_STRING,FOLLOW_2); 
             after(grammarAccess.getElkIdAccess().getIdentifierSTRINGTerminalRuleCall_2_0()); 

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
    // $ANTLR end "rule__ElkId__IdentifierAssignment_2"


    // $ANTLR start "rule__ElkNodeChildren__ChildrenAssignment_1_0"
    // InternalElkGraphJson.g:6288:1: rule__ElkNodeChildren__ChildrenAssignment_1_0 : ( ruleElkNode ) ;
    public final void rule__ElkNodeChildren__ChildrenAssignment_1_0() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:6292:1: ( ( ruleElkNode ) )
            // InternalElkGraphJson.g:6293:2: ( ruleElkNode )
            {
            // InternalElkGraphJson.g:6293:2: ( ruleElkNode )
            // InternalElkGraphJson.g:6294:3: ruleElkNode
            {
             before(grammarAccess.getElkNodeChildrenAccess().getChildrenElkNodeParserRuleCall_1_0_0()); 
            pushFollow(FOLLOW_2);
            ruleElkNode();

            state._fsp--;

             after(grammarAccess.getElkNodeChildrenAccess().getChildrenElkNodeParserRuleCall_1_0_0()); 

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
    // $ANTLR end "rule__ElkNodeChildren__ChildrenAssignment_1_0"


    // $ANTLR start "rule__ElkNodeChildren__ChildrenAssignment_1_1_1"
    // InternalElkGraphJson.g:6303:1: rule__ElkNodeChildren__ChildrenAssignment_1_1_1 : ( ruleElkNode ) ;
    public final void rule__ElkNodeChildren__ChildrenAssignment_1_1_1() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:6307:1: ( ( ruleElkNode ) )
            // InternalElkGraphJson.g:6308:2: ( ruleElkNode )
            {
            // InternalElkGraphJson.g:6308:2: ( ruleElkNode )
            // InternalElkGraphJson.g:6309:3: ruleElkNode
            {
             before(grammarAccess.getElkNodeChildrenAccess().getChildrenElkNodeParserRuleCall_1_1_1_0()); 
            pushFollow(FOLLOW_2);
            ruleElkNode();

            state._fsp--;

             after(grammarAccess.getElkNodeChildrenAccess().getChildrenElkNodeParserRuleCall_1_1_1_0()); 

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
    // $ANTLR end "rule__ElkNodeChildren__ChildrenAssignment_1_1_1"


    // $ANTLR start "rule__ElkNodePorts__PortsAssignment_1_0"
    // InternalElkGraphJson.g:6318:1: rule__ElkNodePorts__PortsAssignment_1_0 : ( ruleElkPort ) ;
    public final void rule__ElkNodePorts__PortsAssignment_1_0() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:6322:1: ( ( ruleElkPort ) )
            // InternalElkGraphJson.g:6323:2: ( ruleElkPort )
            {
            // InternalElkGraphJson.g:6323:2: ( ruleElkPort )
            // InternalElkGraphJson.g:6324:3: ruleElkPort
            {
             before(grammarAccess.getElkNodePortsAccess().getPortsElkPortParserRuleCall_1_0_0()); 
            pushFollow(FOLLOW_2);
            ruleElkPort();

            state._fsp--;

             after(grammarAccess.getElkNodePortsAccess().getPortsElkPortParserRuleCall_1_0_0()); 

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
    // $ANTLR end "rule__ElkNodePorts__PortsAssignment_1_0"


    // $ANTLR start "rule__ElkNodePorts__PortsAssignment_1_1_1"
    // InternalElkGraphJson.g:6333:1: rule__ElkNodePorts__PortsAssignment_1_1_1 : ( ruleElkPort ) ;
    public final void rule__ElkNodePorts__PortsAssignment_1_1_1() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:6337:1: ( ( ruleElkPort ) )
            // InternalElkGraphJson.g:6338:2: ( ruleElkPort )
            {
            // InternalElkGraphJson.g:6338:2: ( ruleElkPort )
            // InternalElkGraphJson.g:6339:3: ruleElkPort
            {
             before(grammarAccess.getElkNodePortsAccess().getPortsElkPortParserRuleCall_1_1_1_0()); 
            pushFollow(FOLLOW_2);
            ruleElkPort();

            state._fsp--;

             after(grammarAccess.getElkNodePortsAccess().getPortsElkPortParserRuleCall_1_1_1_0()); 

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
    // $ANTLR end "rule__ElkNodePorts__PortsAssignment_1_1_1"


    // $ANTLR start "rule__ElkNodeEdges__ContainedEdgesAssignment_1_0"
    // InternalElkGraphJson.g:6348:1: rule__ElkNodeEdges__ContainedEdgesAssignment_1_0 : ( ruleElkEdge ) ;
    public final void rule__ElkNodeEdges__ContainedEdgesAssignment_1_0() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:6352:1: ( ( ruleElkEdge ) )
            // InternalElkGraphJson.g:6353:2: ( ruleElkEdge )
            {
            // InternalElkGraphJson.g:6353:2: ( ruleElkEdge )
            // InternalElkGraphJson.g:6354:3: ruleElkEdge
            {
             before(grammarAccess.getElkNodeEdgesAccess().getContainedEdgesElkEdgeParserRuleCall_1_0_0()); 
            pushFollow(FOLLOW_2);
            ruleElkEdge();

            state._fsp--;

             after(grammarAccess.getElkNodeEdgesAccess().getContainedEdgesElkEdgeParserRuleCall_1_0_0()); 

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
    // $ANTLR end "rule__ElkNodeEdges__ContainedEdgesAssignment_1_0"


    // $ANTLR start "rule__ElkNodeEdges__ContainedEdgesAssignment_1_1_1"
    // InternalElkGraphJson.g:6363:1: rule__ElkNodeEdges__ContainedEdgesAssignment_1_1_1 : ( ruleElkEdge ) ;
    public final void rule__ElkNodeEdges__ContainedEdgesAssignment_1_1_1() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:6367:1: ( ( ruleElkEdge ) )
            // InternalElkGraphJson.g:6368:2: ( ruleElkEdge )
            {
            // InternalElkGraphJson.g:6368:2: ( ruleElkEdge )
            // InternalElkGraphJson.g:6369:3: ruleElkEdge
            {
             before(grammarAccess.getElkNodeEdgesAccess().getContainedEdgesElkEdgeParserRuleCall_1_1_1_0()); 
            pushFollow(FOLLOW_2);
            ruleElkEdge();

            state._fsp--;

             after(grammarAccess.getElkNodeEdgesAccess().getContainedEdgesElkEdgeParserRuleCall_1_1_1_0()); 

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
    // $ANTLR end "rule__ElkNodeEdges__ContainedEdgesAssignment_1_1_1"


    // $ANTLR start "rule__ElkGraphElementLabels__LabelsAssignment_1_0"
    // InternalElkGraphJson.g:6378:1: rule__ElkGraphElementLabels__LabelsAssignment_1_0 : ( ruleElkLabel ) ;
    public final void rule__ElkGraphElementLabels__LabelsAssignment_1_0() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:6382:1: ( ( ruleElkLabel ) )
            // InternalElkGraphJson.g:6383:2: ( ruleElkLabel )
            {
            // InternalElkGraphJson.g:6383:2: ( ruleElkLabel )
            // InternalElkGraphJson.g:6384:3: ruleElkLabel
            {
             before(grammarAccess.getElkGraphElementLabelsAccess().getLabelsElkLabelParserRuleCall_1_0_0()); 
            pushFollow(FOLLOW_2);
            ruleElkLabel();

            state._fsp--;

             after(grammarAccess.getElkGraphElementLabelsAccess().getLabelsElkLabelParserRuleCall_1_0_0()); 

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
    // $ANTLR end "rule__ElkGraphElementLabels__LabelsAssignment_1_0"


    // $ANTLR start "rule__ElkGraphElementLabels__LabelsAssignment_1_1_1"
    // InternalElkGraphJson.g:6393:1: rule__ElkGraphElementLabels__LabelsAssignment_1_1_1 : ( ruleElkLabel ) ;
    public final void rule__ElkGraphElementLabels__LabelsAssignment_1_1_1() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:6397:1: ( ( ruleElkLabel ) )
            // InternalElkGraphJson.g:6398:2: ( ruleElkLabel )
            {
            // InternalElkGraphJson.g:6398:2: ( ruleElkLabel )
            // InternalElkGraphJson.g:6399:3: ruleElkLabel
            {
             before(grammarAccess.getElkGraphElementLabelsAccess().getLabelsElkLabelParserRuleCall_1_1_1_0()); 
            pushFollow(FOLLOW_2);
            ruleElkLabel();

            state._fsp--;

             after(grammarAccess.getElkGraphElementLabelsAccess().getLabelsElkLabelParserRuleCall_1_1_1_0()); 

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
    // $ANTLR end "rule__ElkGraphElementLabels__LabelsAssignment_1_1_1"


    // $ANTLR start "rule__ElkGraphElementProperties__PropertiesAssignment_1_0"
    // InternalElkGraphJson.g:6408:1: rule__ElkGraphElementProperties__PropertiesAssignment_1_0 : ( ruleProperty ) ;
    public final void rule__ElkGraphElementProperties__PropertiesAssignment_1_0() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:6412:1: ( ( ruleProperty ) )
            // InternalElkGraphJson.g:6413:2: ( ruleProperty )
            {
            // InternalElkGraphJson.g:6413:2: ( ruleProperty )
            // InternalElkGraphJson.g:6414:3: ruleProperty
            {
             before(grammarAccess.getElkGraphElementPropertiesAccess().getPropertiesPropertyParserRuleCall_1_0_0()); 
            pushFollow(FOLLOW_2);
            ruleProperty();

            state._fsp--;

             after(grammarAccess.getElkGraphElementPropertiesAccess().getPropertiesPropertyParserRuleCall_1_0_0()); 

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
    // $ANTLR end "rule__ElkGraphElementProperties__PropertiesAssignment_1_0"


    // $ANTLR start "rule__ElkGraphElementProperties__PropertiesAssignment_1_1_1"
    // InternalElkGraphJson.g:6423:1: rule__ElkGraphElementProperties__PropertiesAssignment_1_1_1 : ( ruleProperty ) ;
    public final void rule__ElkGraphElementProperties__PropertiesAssignment_1_1_1() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:6427:1: ( ( ruleProperty ) )
            // InternalElkGraphJson.g:6428:2: ( ruleProperty )
            {
            // InternalElkGraphJson.g:6428:2: ( ruleProperty )
            // InternalElkGraphJson.g:6429:3: ruleProperty
            {
             before(grammarAccess.getElkGraphElementPropertiesAccess().getPropertiesPropertyParserRuleCall_1_1_1_0()); 
            pushFollow(FOLLOW_2);
            ruleProperty();

            state._fsp--;

             after(grammarAccess.getElkGraphElementPropertiesAccess().getPropertiesPropertyParserRuleCall_1_1_1_0()); 

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
    // $ANTLR end "rule__ElkGraphElementProperties__PropertiesAssignment_1_1_1"


    // $ANTLR start "rule__ShapeElement__XAssignment_0_2"
    // InternalElkGraphJson.g:6438:1: rule__ShapeElement__XAssignment_0_2 : ( ruleNumber ) ;
    public final void rule__ShapeElement__XAssignment_0_2() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:6442:1: ( ( ruleNumber ) )
            // InternalElkGraphJson.g:6443:2: ( ruleNumber )
            {
            // InternalElkGraphJson.g:6443:2: ( ruleNumber )
            // InternalElkGraphJson.g:6444:3: ruleNumber
            {
             before(grammarAccess.getShapeElementAccess().getXNumberParserRuleCall_0_2_0()); 
            pushFollow(FOLLOW_2);
            ruleNumber();

            state._fsp--;

             after(grammarAccess.getShapeElementAccess().getXNumberParserRuleCall_0_2_0()); 

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
    // $ANTLR end "rule__ShapeElement__XAssignment_0_2"


    // $ANTLR start "rule__ShapeElement__YAssignment_1_2"
    // InternalElkGraphJson.g:6453:1: rule__ShapeElement__YAssignment_1_2 : ( ruleNumber ) ;
    public final void rule__ShapeElement__YAssignment_1_2() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:6457:1: ( ( ruleNumber ) )
            // InternalElkGraphJson.g:6458:2: ( ruleNumber )
            {
            // InternalElkGraphJson.g:6458:2: ( ruleNumber )
            // InternalElkGraphJson.g:6459:3: ruleNumber
            {
             before(grammarAccess.getShapeElementAccess().getYNumberParserRuleCall_1_2_0()); 
            pushFollow(FOLLOW_2);
            ruleNumber();

            state._fsp--;

             after(grammarAccess.getShapeElementAccess().getYNumberParserRuleCall_1_2_0()); 

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
    // $ANTLR end "rule__ShapeElement__YAssignment_1_2"


    // $ANTLR start "rule__ShapeElement__WidthAssignment_2_2"
    // InternalElkGraphJson.g:6468:1: rule__ShapeElement__WidthAssignment_2_2 : ( ruleNumber ) ;
    public final void rule__ShapeElement__WidthAssignment_2_2() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:6472:1: ( ( ruleNumber ) )
            // InternalElkGraphJson.g:6473:2: ( ruleNumber )
            {
            // InternalElkGraphJson.g:6473:2: ( ruleNumber )
            // InternalElkGraphJson.g:6474:3: ruleNumber
            {
             before(grammarAccess.getShapeElementAccess().getWidthNumberParserRuleCall_2_2_0()); 
            pushFollow(FOLLOW_2);
            ruleNumber();

            state._fsp--;

             after(grammarAccess.getShapeElementAccess().getWidthNumberParserRuleCall_2_2_0()); 

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
    // $ANTLR end "rule__ShapeElement__WidthAssignment_2_2"


    // $ANTLR start "rule__ShapeElement__HeightAssignment_3_2"
    // InternalElkGraphJson.g:6483:1: rule__ShapeElement__HeightAssignment_3_2 : ( ruleNumber ) ;
    public final void rule__ShapeElement__HeightAssignment_3_2() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:6487:1: ( ( ruleNumber ) )
            // InternalElkGraphJson.g:6488:2: ( ruleNumber )
            {
            // InternalElkGraphJson.g:6488:2: ( ruleNumber )
            // InternalElkGraphJson.g:6489:3: ruleNumber
            {
             before(grammarAccess.getShapeElementAccess().getHeightNumberParserRuleCall_3_2_0()); 
            pushFollow(FOLLOW_2);
            ruleNumber();

            state._fsp--;

             after(grammarAccess.getShapeElementAccess().getHeightNumberParserRuleCall_3_2_0()); 

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
    // $ANTLR end "rule__ShapeElement__HeightAssignment_3_2"


    // $ANTLR start "rule__Property__KeyAssignment_0"
    // InternalElkGraphJson.g:6498:1: rule__Property__KeyAssignment_0 : ( rulePropertyKey ) ;
    public final void rule__Property__KeyAssignment_0() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:6502:1: ( ( rulePropertyKey ) )
            // InternalElkGraphJson.g:6503:2: ( rulePropertyKey )
            {
            // InternalElkGraphJson.g:6503:2: ( rulePropertyKey )
            // InternalElkGraphJson.g:6504:3: rulePropertyKey
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
    // InternalElkGraphJson.g:6513:1: rule__Property__ValueAssignment_2_0 : ( ruleStringValue ) ;
    public final void rule__Property__ValueAssignment_2_0() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:6517:1: ( ( ruleStringValue ) )
            // InternalElkGraphJson.g:6518:2: ( ruleStringValue )
            {
            // InternalElkGraphJson.g:6518:2: ( ruleStringValue )
            // InternalElkGraphJson.g:6519:3: ruleStringValue
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
    // InternalElkGraphJson.g:6528:1: rule__Property__ValueAssignment_2_1 : ( ruleNumberValue ) ;
    public final void rule__Property__ValueAssignment_2_1() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:6532:1: ( ( ruleNumberValue ) )
            // InternalElkGraphJson.g:6533:2: ( ruleNumberValue )
            {
            // InternalElkGraphJson.g:6533:2: ( ruleNumberValue )
            // InternalElkGraphJson.g:6534:3: ruleNumberValue
            {
             before(grammarAccess.getPropertyAccess().getValueNumberValueParserRuleCall_2_1_0()); 
            pushFollow(FOLLOW_2);
            ruleNumberValue();

            state._fsp--;

             after(grammarAccess.getPropertyAccess().getValueNumberValueParserRuleCall_2_1_0()); 

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
    // InternalElkGraphJson.g:6543:1: rule__Property__ValueAssignment_2_2 : ( ruleBooleanValue ) ;
    public final void rule__Property__ValueAssignment_2_2() throws RecognitionException {

        		int stackSize = keepStackSize();
        	
        try {
            // InternalElkGraphJson.g:6547:1: ( ( ruleBooleanValue ) )
            // InternalElkGraphJson.g:6548:2: ( ruleBooleanValue )
            {
            // InternalElkGraphJson.g:6548:2: ( ruleBooleanValue )
            // InternalElkGraphJson.g:6549:3: ruleBooleanValue
            {
             before(grammarAccess.getPropertyAccess().getValueBooleanValueParserRuleCall_2_2_0()); 
            pushFollow(FOLLOW_2);
            ruleBooleanValue();

            state._fsp--;

             after(grammarAccess.getPropertyAccess().getValueBooleanValueParserRuleCall_2_2_0()); 

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

    // Delegated rules


 

    public static final BitSet FOLLOW_1 = new BitSet(new long[]{0x0000000000000000L});
    public static final BitSet FOLLOW_2 = new BitSet(new long[]{0x0000000000000002L});
    public static final BitSet FOLLOW_3 = new BitSet(new long[]{0x0400000000000000L});
    public static final BitSet FOLLOW_4 = new BitSet(new long[]{0x1801FFFFFFFF0030L});
    public static final BitSet FOLLOW_5 = new BitSet(new long[]{0x0800000000000000L});
    public static final BitSet FOLLOW_6 = new BitSet(new long[]{0x0800000000000002L});
    public static final BitSet FOLLOW_7 = new BitSet(new long[]{0x0001FFFFFFFF0030L});
    public static final BitSet FOLLOW_8 = new BitSet(new long[]{0x2000000000000000L});
    public static final BitSet FOLLOW_9 = new BitSet(new long[]{0x4000000000000000L});
    public static final BitSet FOLLOW_10 = new BitSet(new long[]{0x1800000000000000L});
    public static final BitSet FOLLOW_11 = new BitSet(new long[]{0x0381FFFFFFFF0030L});
    public static final BitSet FOLLOW_12 = new BitSet(new long[]{0x0000000000000010L});
    public static final BitSet FOLLOW_13 = new BitSet(new long[]{0x007FFFFFFFFF0030L});
    public static final BitSet FOLLOW_14 = new BitSet(new long[]{0x8800000000000010L});
    public static final BitSet FOLLOW_15 = new BitSet(new long[]{0x8C00000000000000L});
    public static final BitSet FOLLOW_16 = new BitSet(new long[]{0x1800000000000030L});
    public static final BitSet FOLLOW_17 = new BitSet(new long[]{0x0000000000000030L});
    public static final BitSet FOLLOW_18 = new BitSet(new long[]{0x00000000000000C0L});
    public static final BitSet FOLLOW_19 = new BitSet(new long[]{0x00000000000060D0L});
    public static final BitSet FOLLOW_20 = new BitSet(new long[]{0xCC0000000000E0D0L});
    public static final BitSet FOLLOW_21 = new BitSet(new long[]{0x440000000000E0D0L});

}