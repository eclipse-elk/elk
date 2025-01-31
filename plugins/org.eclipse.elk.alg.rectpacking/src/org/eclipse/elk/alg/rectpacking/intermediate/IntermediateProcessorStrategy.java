/*******************************************************************************
 * Copyright (c) 2022 sdo and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 * 
 * SPDX-License-Identifier: EPL-2.0 
 *******************************************************************************/
package org.eclipse.elk.alg.rectpacking.intermediate;

import org.eclipse.elk.core.alg.ILayoutProcessor;
import org.eclipse.elk.core.alg.ILayoutProcessorFactory;
import org.eclipse.elk.graph.ElkNode;


/**
 * Definition of available intermediate layout processors for the layered layouter. This enumeration also serves as a
 * factory for intermediate layout processors.
 */
public enum IntermediateProcessorStrategy implements ILayoutProcessorFactory<ElkNode> {

    /*
     * In this enumeration, intermediate layout processors are listed by the earliest slot in which
     * they can sensibly be used. The order in which they are listed is determined by the
     * dependencies on other processors.
     */
    
    /*
     * Before Phase 1
     */
    NODE_SIZE_REORDERER,
    INTERACTIVE_NODE_REORDERER,
    MIN_SIZE_PRE_PROCESSOR,
    
    /*
     * Before Phase 2
     */
    MIN_SIZE_POST_PROCESSOR;

    /* (non-Javadoc)
     * @see org.eclipse.elk.core.alg.ILayoutProcessorFactory#create()
     */
    @Override
    public ILayoutProcessor<ElkNode> create() {
        switch (this) {
        case NODE_SIZE_REORDERER:
            return new NodeSizeReorderer();
        case INTERACTIVE_NODE_REORDERER:
            return new InteractiveNodeReorderer();
        case MIN_SIZE_PRE_PROCESSOR:
            return new MinSizePreProcessor();
        case MIN_SIZE_POST_PROCESSOR:
            return new MinSizePostProcessor();

        default:
            break;
        }
        return null;
    }
    
}