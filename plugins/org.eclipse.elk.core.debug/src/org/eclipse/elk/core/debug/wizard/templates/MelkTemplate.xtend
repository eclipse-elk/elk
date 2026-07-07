/*******************************************************************************
 * Copyright (c) 2019 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.core.debug.wizard.templates

class MelkTemplate {
   static def String buildFileContent(String projectPackage,
       String layoutProviderPath,
       String layoutProvider,
       String idPrefix,
       String algorithmName,
       String label) '''
   package «projectPackage»
   
   import «layoutProviderPath».«layoutProvider»
   import org.eclipse.elk.core.math.ElkPadding
   
   bundle {
       metadataClass «algorithmName»MetadataProvider
       idPrefix «idPrefix»
   }
   
   option reverseInput : boolean {
       label "Reverse Input"
       description
           "True if nodes should be placed in reverse order of their
           appearance in the graph."
       default = false
       targets parents
   }
   
   algorithm «algorithmName»(«layoutProvider») {
       label "«label»"
       description "Please insert a short but informative description here"
       metadataClass options.«algorithmName»Options
       supports reverseInput
       supports org.eclipse.elk.padding = new ElkPadding(10)
       supports org.eclipse.elk.spacing.edgeEdge = 5
       supports org.eclipse.elk.spacing.edgeNode = 10
       supports org.eclipse.elk.spacing.nodeNode = 10
   }'''
}