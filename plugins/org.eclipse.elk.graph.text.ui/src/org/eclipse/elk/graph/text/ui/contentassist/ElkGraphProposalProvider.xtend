/*******************************************************************************
 * Copyright (c) 2016 TypeFox GmbH (http://www.typefox.io) and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.graph.text.ui.contentassist

import com.google.inject.Inject
import com.google.inject.Provider
import java.util.ArrayList
import org.eclipse.elk.core.data.LayoutAlgorithmData
import org.eclipse.elk.core.data.LayoutOptionData
import org.eclipse.elk.graph.text.services.ElkGraphGrammarAccess
import org.eclipse.emf.ecore.EObject
import org.eclipse.jface.viewers.StyledString
import org.eclipse.swt.graphics.Image
import org.eclipse.xtext.Assignment
import org.eclipse.xtext.Keyword
import org.eclipse.xtext.RuleCall
import org.eclipse.xtext.ide.editor.contentassist.ContentAssistContext.Builder
import org.eclipse.xtext.ide.editor.contentassist.ContentAssistEntry
import org.eclipse.xtext.ide.editor.contentassist.IIdeContentProposalAcceptor
import org.eclipse.xtext.ide.editor.contentassist.IdeContentProposalProvider
import org.eclipse.xtext.resource.IEObjectDescription
import org.eclipse.xtext.ui.IImageHelper
import org.eclipse.xtext.ui.editor.contentassist.AbstractContentProposalProvider.NullSafeCompletionProposalAcceptor
import org.eclipse.xtext.ui.editor.contentassist.AbstractJavaBasedContentProposalProvider
import org.eclipse.xtext.ui.editor.contentassist.ContentAssistContext
import org.eclipse.xtext.ui.editor.contentassist.ICompletionProposalAcceptor
import org.eclipse.xtext.util.TextRegion

/**
 * Proposal provider that delegates to the generic IDE implementation.
 */
class ElkGraphProposalProvider extends AbstractJavaBasedContentProposalProvider {
    
    static val MAX_ENTRIES = 1000

    @Inject IdeContentProposalProvider ideProvider
    
    @Inject Provider<Builder> builderProvider
    
    @Inject ElkGraphGrammarAccess grammar
    
    @Inject IImageHelper imageHelper

    override createProposals(ContentAssistContext context, ICompletionProposalAcceptor acceptor) {
        val entries = new ArrayList<Pair<ContentAssistEntry, Integer>>
        val ideAcceptor = new IIdeContentProposalAcceptor {
            override accept(ContentAssistEntry entry, int priority) {
                if (entry !== null)
                    entries += entry -> priority
            }
            override canAcceptMoreProposals() {
                entries.size < MAX_ENTRIES
            }
        }
        ideProvider.createProposals(#[context.getIdeContext], ideAcceptor)
        val uiAcceptor = new NullSafeCompletionProposalAcceptor(acceptor)

        entries.forEach [ p, index |
            val entry = p.key
            val proposal = doCreateProposal(entry.proposal, entry.displayString, entry.image, p.value, context)
            uiAcceptor.accept(proposal)
        ]
    }
    
    override completeAssignment(Assignment object, ContentAssistContext context, ICompletionProposalAcceptor acceptor) {
    }
    
    override completeKeyword(Keyword object, ContentAssistContext context, ICompletionProposalAcceptor acceptor) {
    }
    
    override completeRuleCall(RuleCall object, ContentAssistContext context, ICompletionProposalAcceptor acceptor) {
    }
    
    protected def StyledString getDisplayString(ContentAssistEntry entry) {
        val result = new StyledString(entry.label ?: entry.proposal)
        if (!entry.description.nullOrEmpty)
            result.append(new StyledString(' \u2013 ' + entry.description, StyledString.QUALIFIER_STYLER))
        return result
    }

    private def org.eclipse.xtext.ide.editor.contentassist.ContentAssistContext getIdeContext(ContentAssistContext c) {
        val builder = builderProvider.get
        val replaceRegion = c.replaceRegion
        builder
            .setPrefix(c.prefix)
            .setSelectedText(c.selectedText)
            .setRootModel(c.rootModel)
            .setRootNode(c.rootNode)
            .setCurrentModel(c.currentModel)
            .setPreviousModel(c.previousModel)
            .setCurrentNode(c.currentNode)
            .setLastCompleteNode(c.lastCompleteNode)
            .setOffset(c.offset)
            .setReplaceRegion(new TextRegion(replaceRegion.offset, replaceRegion.length))
            .setResource(c.resource)
        for (grammarElement : c.firstSetGrammarElements) {
            builder.accept(grammarElement)
        }
        return builder.toContext()
    }

    protected def Image getImage(ContentAssistEntry entry) {
        switch source: entry.source {
            IEObjectDescription:
                getImage(source)
            EObject:
                getImage(source)
            LayoutOptionData:
                getImage(source, entry.proposal)
            LayoutAlgorithmData:
                imageHelper.getImage('prop_text.gif')
        }
    }
    
    override protected getImage(EObject eObject) {
        if (eObject instanceof Keyword) {
            val key = switch eObject {
                case grammar.rootNodeAccess.graphKeyword_1_0: 'elkgraph'
                case grammar.elkNodeAccess.nodeKeyword_0: 'elknode'
                case grammar.elkEdgeAccess.edgeKeyword_0: 'elkedge'
                case grammar.elkPortAccess.portKeyword_0: 'elkport'
                case grammar.elkLabelAccess.labelKeyword_0: 'elklabel'
            }
            if (key !== null)
                return imageHelper.getImage(key + '.gif')
        }
        return super.getImage(eObject)
    }
    
    private def getImage(LayoutOptionData option, String value) {
        val key = switch option.type {
            case BOOLEAN:
                if (value == 'false') 'prop_false'
                else 'prop_true'
            case INT: 'prop_int'
            case DOUBLE: 'prop_double'
            case ENUM, case ENUMSET: 'prop_choice'
            default: 'prop_text'
        }
        return imageHelper.getImage(key + '.gif')
    }
    
}
