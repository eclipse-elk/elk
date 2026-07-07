/*******************************************************************************
 * Copyright (c) 2020 Kiel University and others.
 *
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.graph.json.text.ui.contentassist

import com.google.inject.Inject
import com.google.inject.Provider
import java.util.ArrayList
import org.eclipse.jface.viewers.StyledString
import org.eclipse.xtext.ide.editor.contentassist.ContentAssistContext.Builder
import org.eclipse.xtext.ide.editor.contentassist.ContentAssistEntry
import org.eclipse.xtext.ide.editor.contentassist.IIdeContentProposalAcceptor
import org.eclipse.xtext.ide.editor.contentassist.IdeContentProposalProvider
import org.eclipse.xtext.ui.editor.contentassist.AbstractContentProposalProvider.NullSafeCompletionProposalAcceptor
import org.eclipse.xtext.ui.editor.contentassist.ContentAssistContext
import org.eclipse.xtext.ui.editor.contentassist.ICompletionProposalAcceptor
import org.eclipse.xtext.util.TextRegion

/**
 * See https://www.eclipse.org/Xtext/documentation/310_eclipse_support.html#content-assist
 * on how to customize the content assistant.
 */
class ElkGraphJsonProposalProvider extends AbstractElkGraphJsonProposalProvider {
    
    static val MAX_ENTRIES = 1000
    
    @Inject IdeContentProposalProvider ideProvider
    
    @Inject Provider<Builder> builderProvider
    
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
            val proposal = doCreateProposal(entry.proposal, entry.displayString, /*entry.image*/ null, p.value, context)
            uiAcceptor.accept(proposal)
        ]
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
    
    protected def StyledString getDisplayString(ContentAssistEntry entry) {
        val result = new StyledString(entry.label ?: entry.proposal)
        if (!entry.description.nullOrEmpty)
            result.append(new StyledString(' \u2013 ' + entry.description, StyledString.QUALIFIER_STYLER))
        return result
    }
    
}
